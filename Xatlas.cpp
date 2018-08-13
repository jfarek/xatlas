#include "Xatlas.hpp"
#include "Bam.hpp"
#include "CoverageCounter.hpp"
#include "EventScanner.hpp"
#include "IndelEvent.hpp"
#include "Logit.hpp"
#include "ReferenceSequence.hpp"
#include "SnpEvent.hpp"
#include "VcfWriter.hpp"
#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include <algorithm>
#include <cstring>
#include <fstream>
#include <getopt.h>
#include <iostream>
#include <sstream>

/* ------===< xAtlas >===------ */

/**
 * Yeah the code is a mess but you should see what I was given to start with ...
 *
 * TODO clean up the code
 * TODO fix cutoffs/thresholds/etc
 * TODO deduplicate the process_*() functions
 * TODO untangle snp and indel processing threads from read buffer thread so we can use semaphores instead of phtread barriers
 * TODO options for calling only snps/only indels
 * TODO more configuration options in general
 */

static const char *help = "\
Required arguments:\n\
    -r, --ref REF           Reference genome in FASTA format\n\
    -i, --in IN             Sorted and indexed input BAM or CRAM file\n\
    -s, --sample-name SN    Sample name to use in the output VCF file\n\
    -p, --prefix PFX        Output VCF file prefix\n\
\n\
Options:\n\
    -P, --multithread               Read alignment file and process records in separate threads\n\
    -t, --num-hts-threads NTHREAD   Number of HTSlib decompression threads to spawn\n\
    -c, --capture-bed BED           BED file of regions to process\n\
    -v, --min-p-value               Minimum logit P-value to report variants\n\
    -m, --min-snp-mapq MAPQ         Minimum read mapping quality for calling SNPs\n\
    -n, --min-indel-mapq MAPQ       Minimum read mapping quality for calling indels\n\
    -M, --max-coverage COV          Maximum coverage for calling variants normally\n\
    -A, --block-abs-lim LIM         gVCF non-variant block absolute range limit\n\
    -R, --block-rel-lim LIM         gVCF non-variant block relative range limit coefficient\n\
    -g, --gvcf                      Include non-variant gVCF blocks in output VCF file\n\
    -z, --bgzf                      Write output in bgzip-compressed VCF format\n\
    -S, --snp-logit-params FILE     File with intercept and coefficients for SNP logit model\n\
    -I, --indel-logit-params FILE   File with intercept and coefficients for indel logit model\n\
    -F, --enable-strand-filter      Enable SNP filter for single-strandedness\n\
    -h, --help                      Show this help\n\
";

struct args {
    opts_s *opts;

    bam1_t **rec_buff;
    bam1_t **read_start;
    bam1_t **process_start;
    bam1_t **curr_rec;

    uint32_t read_section;
    uint32_t buff_size;
    uint32_t num_sections;
    uint32_t section_size;
    uint32_t process_size;
    uint32_t num_pieces;
    uint32_t region_idx;

    Bam *bam;
    hts_itr_t *iter;

    EventScanner *events;
    ReferenceSequence *refseq;
    CoverageCounter *coverages;
    VcfWriter *writer;
    bed_list_t *bedlist;
    qual_t min_indel_mapq;
    qual_t min_snp_mapq;
    bool more_chrs;

    cigar_list_t *cigar_list_read;
    cigar_list_t *cigar_list_process;

    bed_coord_list_t *regions;
    size_t idx;

    void (*regions_itr_func)(struct args *);

#ifdef USE_PTHREAD
    pthread_barrier_t barriers[4];
#endif /* USE_PTHREAD */
};
typedef struct args args_s;

void clear_processed_cigars(args_s *args)
{
    // clear cigar list buffer up to last snp or indel call
    int32_t last_call = std::min(args->events->_last_call_snp, args->events->_last_call_indel);
    auto cigar_it0 = (*args->cigar_list_process).begin();
    auto cigar_it = (*args->cigar_list_process).begin();
    auto cigar_end = (*args->cigar_list_process).end();

    while (cigar_it < cigar_end && cigar_it->first.second < last_call) {
        ++cigar_it;
    }
    args->cigar_list_process->erase(cigar_it0, cigar_it);

    for (auto &c : *args->cigar_list_read) {
        args->cigar_list_process->push_back(c);
    }
    args->cigar_list_read->clear();
}

uint32_t read_bam_section(args_s *args)
{
    uint32_t num_read = 0;
    uint32_t n_cigar, *cigar, *cigar0;
    const int16_t filtered_flags = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    bam1_t *curr_rec;

    while (num_read < args->section_size &&
           sam_itr_next(args->bam->_sf, args->iter, *args->curr_rec) >= 0)
    {
        curr_rec = *args->curr_rec;
        // filters common to snp and indel
        if ((curr_rec->core.flag & filtered_flags) == 0 &&
            curr_rec->core.n_cigar > 0 &&
            bam_aux2i(bam_aux_get(*args->curr_rec, "NM")) != -1)
        {
            std::vector< uint32_t > cigar_vec;
            n_cigar = curr_rec->core.n_cigar;
            cigar0 = bam_get_cigar(curr_rec);
            cigar = cigar0;

            for (uint32_t i = 0; i < n_cigar; ++i) {
                cigar_vec.push_back(*cigar);
                ++cigar;
            }

            // TODO generalize this cigar list buffer thing so that it acts as a bam1_t buffer
            // so we can maybe ditch CovCounter in scanning for variants and recording coverage
            args->cigar_list_read->push_back(
                std::make_pair(
                    std::make_pair(
                        curr_rec->core.pos,
                        curr_rec->core.pos + bam_cigar2rlen(n_cigar, cigar0)),
                    cigar_vec));

            ++args->curr_rec;
            ++num_read;
        }
    }

    return num_read;
}

void process_records_snp_section(args_s *args, bed_coord_t &region)
{
    uint32_t num_read = 0;
    bam1_t *curr_rec;
    bam1_t **next_rec = args->process_start;

    while (num_read < args->process_size) {
        curr_rec = *next_rec;
        ++next_rec;
        ++num_read;

        if (curr_rec->core.qual >= args->min_snp_mapq &&
            !args->coverages->test_high_coverage(curr_rec->core.pos, COVSIDX_SNP_DP, curr_rec->core.l_qseq))
        {
            args->events->collect_snps(curr_rec, *args->refseq, *args->coverages, region);
            if (!args->events->_snps.empty()) {
                args->writer->print_snp_buffer(curr_rec->core.pos, region);
            }
        }
    }
}

void process_records_indel_section(args_s *args, bed_coord_t &region)
{
    uint32_t num_read = 0;
    bam1_t *curr_rec;
    bam1_t **next_rec = args->process_start;

    while (num_read < args->process_size) {
        curr_rec = *next_rec;
        ++next_rec;
        ++num_read;

        if (curr_rec->core.qual >= args->min_indel_mapq &&
            !args->coverages->test_high_coverage(curr_rec->core.pos, COVSIDX_INDEL_DP, curr_rec->core.l_qseq))
        {
            args->events->collect_indels(curr_rec, *args->refseq, *args->coverages, region);
            if (!args->events->_indels.empty()) {
                args->writer->print_indel_buffer(curr_rec->core.pos, region);
            }
        }
    }
}

#ifdef USE_PTHREAD
void *process_records_indel(void *arg)
{
    args_s *args = (args_s *)arg;

    // for each chrom
    while (args->more_chrs) {
        pthread_barrier_wait(&args->barriers[0]);

        bed_coord_t &region = (*args->regions)[args->idx];

        goto skip_first_pass_indels;
        do {
            process_records_indel_section(args, region);

        skip_first_pass_indels:
            pthread_barrier_wait(&args->barriers[1]);
            pthread_barrier_wait(&args->barriers[2]);

        } while (args->process_size > 0);

        // final variants and gvcf block
        args->writer->print_indel_buffer(args->opts->last_call_max, region);

        if (args->opts->gvcfoutput) {
            if (args->events->_last_call_indel < region.first) {
                args->events->_last_call_indel = region.first - 1;
            }
            args->writer->print_gvcf_span(region.second + 1,
                args->events->_last_call_indel,
                region.second,
                false,
                args->writer->_prev_del,
                COVSIDX_INDEL_VR,
                COVSIDX_INDEL_RR,
                COVSIDX_INDEL_DP);
        }

        args->writer->_prev_del = false;

        pthread_barrier_wait(&args->barriers[3]);
    }

    pthread_exit(nullptr);
}

void *process_records_snp(void *arg)
{
    args_s *args = (args_s *)arg;

    // for each chrom
    while (args->more_chrs) {
        pthread_barrier_wait(&args->barriers[0]);

        bed_coord_t &region = (*args->regions)[args->idx];

        goto skip_first_pass_snps;
        do {
            process_records_snp_section(args, region);

        skip_first_pass_snps:
            pthread_barrier_wait(&args->barriers[1]);
            pthread_barrier_wait(&args->barriers[2]);

        } while (args->process_size > 0);

        // final variants and gvcf block
        args->writer->print_snp_buffer(args->opts->last_call_max, region);

        if (args->opts->gvcfoutput) {
            if (args->events->_last_call_snp < region.first) {
                args->events->_last_call_snp = region.first - 1;
            }
            args->writer->print_gvcf_span(region.second + 1,
                args->events->_last_call_snp,
                region.second,
                true,
                false,
                COVSIDX_SNP_VR,
                COVSIDX_SNP_RR,
                COVSIDX_SNP_DP);
        }

        pthread_barrier_wait(&args->barriers[3]);
    }

    pthread_exit(nullptr);
}

void read_bam(args_s *args)
{
    uint32_t num_read;

    pthread_barrier_wait(&args->barriers[0]);

    do {
        args->curr_rec = args->read_start;
        num_read = read_bam_section(args);

        pthread_barrier_wait(&args->barriers[1]);

        args->process_start = args->read_start;
        args->process_size = num_read;

        ++args->read_section;
        args->read_section %= args->num_sections;

        args->read_start = (args->read_section == 0)
                               ? args->rec_buff
                               : args->curr_rec;

        clear_processed_cigars(args);

        pthread_barrier_wait(&args->barriers[2]);
    } while (args->process_size > 0);

    ++args->region_idx;
    args->more_chrs = (args->region_idx < args->num_pieces);

    pthread_barrier_wait(&args->barriers[3]);
}
#endif /* USE_PTHREAD */

/**
 * Single-threaded reading and processing
 */
void read_and_process_bam(args_s *args)
{
    bed_coord_t &region = (*args->regions)[args->idx];

    do {
        args->curr_rec = args->read_start;
        args->process_size = read_bam_section(args);

        args->process_start = args->read_start;

        ++args->read_section;
        args->read_section %= args->num_sections;

        args->read_start = (args->read_section == 0)
                               ? args->rec_buff
                               : args->curr_rec;

        clear_processed_cigars(args);

        if (args->process_size > 0) {
            process_records_snp_section(args, region);
            process_records_indel_section(args, region);
        } else {
            break;
        }
    } while (1);

    // snp final variants and gvcf block
    args->writer->print_snp_buffer(args->opts->last_call_max, region);

    if (args->opts->gvcfoutput) {
        if (args->events->_last_call_snp < region.first) {
            args->events->_last_call_snp = region.first - 1;
        }
        args->writer->print_gvcf_span(
            region.second + 1,
            args->events->_last_call_snp,
            region.second,
            true,
            false,
            COVSIDX_SNP_VR,
            COVSIDX_SNP_RR,
            COVSIDX_SNP_DP);
    }

    // indel final variants and gvcf block
    args->writer->print_indel_buffer(args->opts->last_call_max, region);

    if (args->opts->gvcfoutput) {
        if (args->events->_last_call_indel < region.first) {
            args->events->_last_call_indel = region.first - 1;
        }
        args->writer->print_gvcf_span(
            region.second + 1,
            args->events->_last_call_indel,
            region.second,
            false,
            args->writer->_prev_del,
            COVSIDX_INDEL_VR,
            COVSIDX_INDEL_RR,
            COVSIDX_INDEL_DP);
    }

    args->writer->_prev_del = false;
}

void iterate_regions(args_s *args)
{
    bed_coord_list_t regions;
    std::vector< std::string > chrom_regions;
    args->region_idx = 0;

    // for each chrom
    for (auto &chrom : *args->bedlist) {
        const char *reg = chrom.first.c_str();
        regions = chrom.second;

        args->writer->set_region(reg);
        args->refseq->set_region(reg);
        args->events->reset();
        args->coverages->reset();
        args->cigar_list_read->clear();
        args->cigar_list_process->clear();
        chrom_regions.reserve(regions.size());
        chrom_regions.clear();

        for (const auto &region : regions) {
            if (region.second >= region.first) {
                std::stringstream x;
                x << reg << ":" << (region.first + 1) << "-" << (region.second + 1);
                chrom_regions.push_back(x.str());
            }
        }

        args->read_start = args->rec_buff;
        args->process_start = args->rec_buff;
        args->read_section = 0;
        args->process_size = 0;
        args->regions = &regions;

        // for each piece in chrom
        for (size_t idx = 0; idx < chrom_regions.size(); ++idx) {
            std::cerr << "Processing records for region " 
                      << reg << ":" << (regions[idx].first + 1) << "-" << (regions[idx].second + 1) << std::endl;

            args->events->_last_call_indel = regions[idx].first;
            args->events->_last_call_snp = regions[idx].first;
            args->bam->set_iter(chrom_regions[idx].c_str());
            args->iter = args->bam->_iter;
            args->idx = idx;

            args->regions_itr_func(args);

            if (args->events->_snps.size() + args->events->_indels.size() > 0) {
                if (!args->events->_snps.empty()) {
                    std::cerr << "Error: " << args->events->_snps.size() << " remaining unprocessed SNPs." << std::endl;
                }
                if (!args->events->_indels.empty()) {
                    std::cerr << "Error: " << args->events->_indels.size() << " remaining unprocessed INDELs." << std::endl;
                }
                exit(EXIT_FAILURE);
            }
        }
    }
}

#ifdef USE_PTHREAD
void *iterate_regions_pt(void *arg)
{
    iterate_regions((args_s *)arg);

    pthread_exit(nullptr);
}
#endif /* USE_PTHREAD */

void read_logit_params(logit_params_s *logit_params, const char *fn, bool snp)
{
    const size_t max_values = 6;
    std::string buffer;
    std::vector< double > values;

    std::ifstream ifs(fn);
    if (!ifs) {
        std::cerr << "Warning: failed to read logistic regression parameter file \"" << fn << "\", "
                  << "using default values." << std::endl;
        return;
    }

    for (size_t i = 0; i < max_values && std::getline(ifs, buffer); ++i) {
        values.push_back(std::strtof(buffer.c_str(), nullptr));
    }

    if (snp) {
        logit_params->snp->intercept = values[0];
        logit_params->snp->ratio_score = values[1];
        logit_params->snp->base_qual = values[2];
        logit_params->snp->mean_avnqs = values[3];
        logit_params->snp->rel_pos = values[4];
        logit_params->snp->titv = values[5];
    } else {
        logit_params->indel->intercept = values[0];
        logit_params->indel->ratio_score = values[1];
        logit_params->indel->strand_dir = values[2];
        logit_params->indel->mean_avnqs = values[3];
        logit_params->indel->seq_entropy = values[4];
        logit_params->indel->mean_var_rate = values[5];
    }

    ifs.close();
}

void usage()
{
    std::cerr << XATLAS_NAME << " v" << XATLAS_VERSION << std::endl
              << std::endl
              << help << std::endl;
}

int main(int argc, char **argv)
{
    bool do_multithreading = false;
    char *ref = nullptr;
    char *bam_fn = nullptr;
    char *pfx = nullptr;
    char *sample_name = nullptr;
    qual_t min_indel_mapq, min_snp_mapq;
    int tmp_min_indel_mapq = 1;
    int tmp_min_snp_mapq = 1;
    uint8_t num_hts_threads = 1;

    coverage_t max_cov = 16383;
    double snp_max_sub = 0.05;
    double snp_max_indel = 0.05;

    // Options

    opts_s opts;

    opts.gvcfoutput = false;
    opts.bgzf = false;
    opts.dumpsnp = false;
    opts.dumpindel = false;
    //opts.report_read_end_score = false;
    opts.capturebed = nullptr;
    opts.block_abs_lim = 3;
    opts.block_rel_lim = 0.3;
    opts.block_rel_min = 1.0;
    opts.last_call_max = INT32_MAX;
    opts.high_cov_cutoff = 8000;
    opts.min_pr = 0.25;

    opts.snp_min_dp = 6;
    opts.snp_min_vr = 2;
    opts.enable_snp_strand_cutoff = false;
    opts.snp_strand_cutoff = 16;
    opts.snp_strand_ratio_cutoff = 0.01;
    opts.snp_near_end_bases = 3;
    //opts.snp_het_min = 0.25; //0.1;
    //opts.snp_het_max = 0.75; //0.9;

    opts.indel_min_depth = 5;
    opts.indel_min_var_reads = 2;
    opts.indel_min_var_ratio = 0.06;
    opts.indel_max_near_read_end_ratio = 0.8;
    opts.indel_het_min = 0.25; //0.06;
    opts.indel_het_max = 0.75; //0.6;
    //opts.indel_strand_dir_filter = false;

    snp_logit_params_s logit_params_snp;
    logit_params_snp.intercept = -6.66404;
    logit_params_snp.ratio_score = 11.1192;
    logit_params_snp.base_qual = 0.25579;
    logit_params_snp.mean_avnqs = -0.12896;
    logit_params_snp.rel_pos = -0.69106;
    logit_params_snp.titv = 0.4851;

    indel_logit_params_s logit_params_indel;
    logit_params_indel.intercept = -7.1085;
    logit_params_indel.ratio_score = 6.22804;
    logit_params_indel.strand_dir = 2.21407;
    logit_params_indel.mean_avnqs = 0.07777;
    logit_params_indel.seq_entropy = 0.1479;
    logit_params_indel.mean_var_rate = -2.13305;

    logit_params_s logit_params;
    logit_params.snp = &logit_params_snp;
    logit_params.indel = &logit_params_indel;

    // Runtime options

    int c, optidx;
    static struct option long_options[] = {
        {"ref", 1, nullptr, 0},                     // r
        {"in", 1, nullptr, 0},                      // i
        {"sample-name", 1, nullptr, 0},             // s
        {"prefix", 1, nullptr, 0},                  // p
        {"multithread", 0, nullptr, 0},             // P
        {"num-hts-threads", 1, nullptr, 0},         // t
        {"capture-bed", 1, nullptr, 0},             // c
        {"min-p-value", 1, nullptr, 0},             // v
        {"min-snp-mapq", 1, nullptr, 0},            // m
        {"min-indel-mapq", 1, nullptr, 0},          // n
        {"max-coverage", 1, nullptr, 0},            // M
        {"block-abs-lim", 1, nullptr, 0},           // A
        {"block-rel-lim", 1, nullptr, 0},           // R
        {"gvcf", 0, nullptr, 0},                    // g
        {"bgzf", 0, nullptr, 0},                    // z
        {"snp-logit-params", 1, nullptr, 0},        // S
        {"indel-logit-params", 1, nullptr, 0},      // I
        {"enable-strand-filter", 0, nullptr, 0},    // F
        {"help", 0, nullptr, 0},                    // h
        // hidden options
        {"dump-snp", 0, nullptr, 0},       // Z
        {"dump-indel", 0, nullptr, 0},     // Y
        //{"read-end-score", 0, nullptr, 0}, // X
        {nullptr, 0, nullptr, 0}};
    const char *short_options = "0r:i:s:p:Pt:c:v:m:n:M:A:R:gzS:I:FhZY";
    const char *shorter_options = "rispPtcvmnMARgzSIFhZY";

    if (argc == 1) {
        usage();
        return EXIT_SUCCESS;
    }

    while ((c = getopt_long(argc, argv, short_options, long_options, &optidx)) != -1) {
        if (c == 0) {
            c = shorter_options[optidx];
        }

        switch (c) {
        case 'r': // ref
            ref = optarg;
            break;
        case 'i': // in
            bam_fn = optarg;
            break;
        case 's': // sample-name
            sample_name = optarg;
            break;
        case 'p': // prefix
            pfx = optarg;
            break;
        case 'P': // prefix
            do_multithreading = true;
            break;
        case 't': // num-hts-threads
            num_hts_threads = (uint8_t)std::atoi(optarg);
            break;
        case 'c': // capture-bed
            opts.capturebed = optarg;
            break;
        case 'v': // min-p-value
            opts.min_pr = std::strtod(optarg, nullptr);
            break;
        case 'm': // min-snp-mapq
            tmp_min_snp_mapq = (qual_t)std::atoi(optarg);
            break;
        case 'n': // min-indel-mapq
            tmp_min_indel_mapq = (qual_t)std::atoi(optarg);
            break;
        case 'M': // max-coverage
            opts.high_cov_cutoff = std::atoi(optarg);
            break;
        case 'A': // block-abs-lim
            opts.block_abs_lim = std::atoi(optarg);
            break;
        case 'R': // block-rel-lim
            opts.block_rel_lim = std::strtod(optarg, nullptr);
            break;
        case 'g': // gvcf
            opts.gvcfoutput = true;
            break;
        case 'z': // bgzf
            opts.bgzf = true;
            break;
        case 'S': // snp-logit-params
            read_logit_params(&logit_params, optarg, true);
            break;
        case 'I': // indel-logit-params
            read_logit_params(&logit_params, optarg, false);
            break;
        case 'F': // enable-strand-filter
            opts.enable_snp_strand_cutoff = true;
            break;
        case 'Z': // dump-snp
            opts.dumpsnp = true;
            break;
        case 'Y': // dump-indel
            opts.dumpindel = true;
            break;
        /*
        case 'X': // read-end-score
            opts.report_read_end_score = true;
            break;
        */
        case 'h': // help
            usage();
            return EXIT_SUCCESS;
        default:
            return EXIT_FAILURE;
        }
    }

    snprintf(opts.block_label, 24, "BLOCKAVG_min%dp%da", (int)(100 * opts.block_rel_lim), (int)opts.block_abs_lim);

    // Begin

    std::cerr << "Running " << XATLAS_NAME << " v" << XATLAS_VERSION << std::endl
              << "SNP Minimum coverage: " << opts.snp_min_dp << std::endl
              << "SNP/INDEL High coverage: " << opts.high_cov_cutoff << std::endl
              << "SNP/INDEL Maximum coverage: " << max_cov << std::endl
              << "gVCF block absolute limit: " << opts.block_abs_lim << std::endl
              << "gVCF block relative limit: " << opts.block_rel_lim << std::endl
              << "gVCF block label: " << opts.block_label << std::endl;

    if (opts.high_cov_cutoff > 8191) {
        std::cerr << "Maximum coverage is too high" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (bam_fn != nullptr) {
        std::cerr << "Input alignment file: " << bam_fn << std::endl;
    } else {
        std::cerr << "No input alignment file given" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (ref != nullptr) {
        std::cerr << "Reference file: " << ref << std::endl;
    } else {
        std::cerr << "No reference file given" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (sample_name != nullptr) {
        std::cerr << "Sample name: " << sample_name << std::endl;
    } else {
        std::cerr << "No sample name given" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (tmp_min_indel_mapq < 0 || tmp_min_indel_mapq > 60) {
        std::cerr << "Minimum INDEL mapping quality must be in the range [0, 60]" << std::endl;
        exit(EXIT_FAILURE);
    } else if (tmp_min_snp_mapq < 0 || tmp_min_snp_mapq > 60) {
        std::cerr << "Minimum SNP mapping quality must be in the range [0, 60]" << std::endl;
        exit(EXIT_FAILURE);
    } else if (tmp_min_snp_mapq < tmp_min_indel_mapq) {
        std::cerr << "Setting lower minimum mapping qual for SNPs has not been implemented" << std::endl;
        exit(EXIT_FAILURE);
    }

    min_snp_mapq = tmp_min_snp_mapq;
    min_indel_mapq = tmp_min_indel_mapq;

    // Files

    if (pfx == nullptr) {
        std::cerr << "No filename prefix given" << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cerr << "Prefix: " << pfx << std::endl;
    std::cerr << "Minimum SNP read mapping quality: " << (short)min_snp_mapq << std::endl;
    std::cerr << "Minimum INDEL read mapping quality: " << (short)min_indel_mapq << std::endl;

    Bam bam(bam_fn, ref, num_hts_threads);
    switch (bam._status) {
    case BAM_BAD_FILE:
        std::cerr << "Failed to load alignment file \"" << bam_fn << "\"" << std::endl;
        exit(EXIT_FAILURE);
        break;
    case BAM_BAD_HEADER:
        std::cerr << "Failed to load header for \"" << bam_fn << "\"" << std::endl;
        exit(EXIT_FAILURE);
        break;
    case BAM_BAD_INDEX:
        std::cerr << "Failed to load index for \"" << bam_fn << "\"" << std::endl;
        exit(EXIT_FAILURE);
        break;
    case BAM_BAD_REF_INDEX:
        std::cerr << "Failed to load reference index for \"" << ref << "\"" << std::endl;
        exit(EXIT_FAILURE);
        break;
    default:
        break;
    }

    std::cerr << "Found " << bam._hdr->n_targets << " regions" << std::endl;

    std::string pfx_str(pfx);
    std::string sfx_str(opts.bgzf ? "vcf.gz" : "vcf");
    std::string indel_fn(pfx_str + "_indel." + sfx_str), snp_fn(pfx_str + "_snp." + sfx_str);

    vcfFile *indel_fp, *snp_fp;
    const char *write_mode = opts.bgzf ? "wz" : "w";
    if ((indel_fp = hts_open(indel_fn.c_str(), write_mode)) == nullptr) {
        std::cerr << "Unable to open output INDEL VCF file" << std::endl;
        exit(EXIT_FAILURE);
    }
    if ((snp_fp = hts_open(snp_fn.c_str(), write_mode)) == nullptr) {
        std::cerr << "Unable to open output SNP VCF file" << std::endl;
        exit(EXIT_FAILURE);
    }

    bcf_hdr_t *indel_hdr, *snp_hdr;
    if ((indel_hdr = bcf_hdr_init("w")) == nullptr) {
        std::cerr << "Unable to create INDEL header" << std::endl;
        exit(EXIT_FAILURE);
    }
    bcf_hdr_set_version(indel_hdr, "VCFv4.3");
    if ((snp_hdr = bcf_hdr_init("w")) == nullptr) {
        std::cerr << "Unable to create SNP header" << std::endl;
        exit(EXIT_FAILURE);
    }
    bcf_hdr_set_version(snp_hdr, "VCFv4.3");

    // Which regions

    regions_list_t regions;
    regions.reserve(bam._hdr->n_targets);

    for (int i = 0; i < bam._hdr->n_targets; ++i) {
        regions.push_back(std::make_pair(std::string(bam._hdr->target_name[i]), bam._hdr->target_len[i]));
    }

    ReferenceSequence refseq(ref);
    bed_list_t bedlist;
    uint32_t num_pieces = 0;

    if (opts.capturebed != nullptr) {
        bed_coord_map_t *bedmap = new bed_coord_map_t;
        std::string tmp;

        std::cerr << "Using " << opts.capturebed << " for input sequence list" << std::endl;
        std::ifstream ifs(opts.capturebed);
        if (!ifs) {
            std::cerr << "Problem reading " << opts.capturebed << std::endl;
            exit(EXIT_FAILURE);
        }

        char *tname = new char[64];
        int32_t start, end;
        while (std::getline(ifs, tmp)) {
            if ((std::sscanf(tmp.c_str(), "%63s %d %d", tname, &start, &end)) < 3) {
                continue;
            }
            (*bedmap)[std::string(tname)].push_back(std::make_pair(start, end - 1));
            ++num_pieces;
        }
        delete[] tname;

        ifs.close();

        for (const auto &region : regions) {
            if (bedmap->count(region.first) > 0) {
                bedlist.push_back(std::make_pair(region.first, (*bedmap)[region.first]));
            }
        }

        delete bedmap;
    } else {
        for (const auto &region : regions) {
            std::vector< bed_coord_t > vec;
            vec.push_back(std::make_pair(0, faidx_seq_len(refseq._fai, region.first.c_str()) - 1));
            bedlist.push_back(std::make_pair(region.first, vec));
            ++num_pieces;
        }
    }

    if (bedlist.empty()) {
        std::cerr << "No regions to process" << std::endl;
        exit(EXIT_FAILURE);
    }

    // max region length
    int32_t max_len = 0;
    for (const auto &bed_coord : bedlist) {
        for (const auto &coord : bed_coord.second) {
            if (coord.second > max_len) {
                max_len = coord.second;
            }
        }
    }

    // Setup

    args_s args;

    args.opts = &opts;
    args.num_sections = 4;
    args.section_size = 0x4000;
    args.min_indel_mapq = min_indel_mapq;
    args.min_snp_mapq = min_snp_mapq;
    args.buff_size = args.num_sections * args.section_size;
    args.rec_buff = new bam1_t *[args.buff_size];
    for (uint32_t k = 0; k < args.buff_size; ++k) {
        args.rec_buff[k] = bam_init1();
    }
    args.num_pieces = num_pieces;

    cigar_list_t cigar_list_read;
    cigar_list_t cigar_list_process;

    EventScanner events(opts.snp_near_end_bases, snp_max_sub, snp_max_indel);
    CoverageCounter coverages((size_t)max_len, max_cov);
    VcfWriter writer(coverages, refseq, events, cigar_list_process, indel_hdr, snp_hdr, indel_fp, snp_fp, &opts, &logit_params);

    writer.setup_vcf(sample_name, argc, argv, XATLAS_VERSION, ref, regions);
    bcf_hdr_write(indel_fp, indel_hdr);
    bcf_hdr_write(snp_fp, snp_hdr);

    args.bam = &bam;
    args.events = &events;
    args.refseq = &refseq;
    args.coverages = &coverages;
    args.writer = &writer;
    args.bedlist = &bedlist;
    args.cigar_list_read = &cigar_list_read;
    args.cigar_list_process = &cigar_list_process;
    args.more_chrs = true;

    if (do_multithreading) {
#ifdef USE_PTHREAD
        args.regions_itr_func = read_bam;

        for (size_t i = 0; i < 4; ++i) {
            pthread_barrier_init(&args.barriers[i], nullptr, 3);
        }
        pthread_t read_bam_thread, process_records_snp_thread, process_records_indel_thread;
        void *status;

        if (pthread_create(&read_bam_thread, nullptr, iterate_regions_pt, &args) != 0) {
            std::cerr << "Error: Unable to create read_bam_thread" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_create(&process_records_indel_thread, nullptr, process_records_indel, &args) != 0) {
            std::cerr << "Error: Unable to create process_records_indel_thread" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_create(&process_records_snp_thread, nullptr, process_records_snp, &args) != 0) {
            std::cerr << "Error: Unable to create process_records_snp_thread" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_join(read_bam_thread, &status) != 0) {
            std::cerr << "Error: Unable to join read_bam_thread" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_join(process_records_indel_thread, &status) != 0) {
            std::cerr << "Error: Unable to join process_records_indel_thread" << std::endl;
            exit(EXIT_FAILURE);
        }

        if (pthread_join(process_records_snp_thread, &status) != 0) {
            std::cerr << "Error: Unable to join process_records_snp_thread" << std::endl;
            exit(EXIT_FAILURE);
        }

        // Cleanup

        for (size_t i = 0; i < 4; ++i) {
            pthread_barrier_destroy(&args.barriers[i]);
        }
    } else {
#else
        std::cerr << "Multithreading not supported in this build of xAtlas" << std::endl;
    }
#endif /* USE_PTHREAD */

        args.regions_itr_func = read_and_process_bam;
        iterate_regions(&args);

#ifdef USE_PTHREAD
    }
#endif /* USE_PTHREAD */

    for (uint32_t k = 0; k < args.buff_size; ++k) {
        bam_destroy1(args.rec_buff[k]);
    }
    delete[] args.rec_buff;

    hts_close(indel_fp);
    hts_close(snp_fp);
    bcf_hdr_destroy(indel_hdr);
    bcf_hdr_destroy(snp_hdr);

    std::cerr << "Finished" << std::endl;

    return EXIT_SUCCESS;
}
