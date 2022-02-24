#include "xatlas/VcfWriter.h"
#include "xatlas/GvcfBlock.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <sstream>

/**
 * Main VCF writing facilities
 */

VcfWriter::VcfWriter()
{
}

VcfWriter::VcfWriter(
        CoverageCounter &coverages,
        ReferenceSequence &sequences,
        EventScanner &events,
        cigar_list_t &cigar_list,
        bcf_hdr_t *indel_hdr,
        bcf_hdr_t *snp_hdr,
        htsFile *indel_fp,
        htsFile *snp_fp,
        const opts_s *opts,
        const logit_params *logit_params)
    : _coverages(&coverages),
      _sequences(&sequences),
      _events(&events),
      _cigar_list(&cigar_list),
      _region(""),
      _indel_rid(-1),
      _snp_rid(-1),
      _indel_hdr(indel_hdr),
      _snp_hdr(snp_hdr),
      _indel_fp(indel_fp),
      _snp_fp(snp_fp),
      _opts(opts),
      _logit_params(logit_params),
      _prev_del_end(-1),
      _prev_del(false)
{
    if (((_indel_rec = bcf_init1()) != nullptr) &&
        ((_indel_gvcf_rec = bcf_init1()) != nullptr) &&
        ((_snp_rec = bcf_init1()) != nullptr) &&
        ((_snp_gvcf_rec = bcf_init1()) != nullptr))
    {
        _gvcf_label = opts->block_label;
        _indel_pass_int = bcf_hdr_id2int(indel_hdr, BCF_DT_ID, "PASS");
        _snp_pass_int = bcf_hdr_id2int(snp_hdr, BCF_DT_ID, "PASS");
        _prev_del_gt[0] = -1;
        _prev_del_gt[1] = -1;
        _gt_missing = bcf_int32_missing;

        std::strncpy(_snp_alleles, ".,.", 4);
        std::strncpy(_snp_alleles_gvcf, ".,.", 4);
        std::strncpy(_indel_alleles_gvcf, ".,.", 4);
    }
}

VcfWriter::VcfWriter(const VcfWriter &) = default;

VcfWriter::~VcfWriter()
{
    bcf_destroy1(_indel_rec);
    bcf_destroy1(_indel_gvcf_rec);
    bcf_destroy1(_snp_rec);
    bcf_destroy1(_snp_gvcf_rec);
}

void VcfWriter::set_region(const char *reg)
{
    _region = reg;
    _indel_rid = bcf_hdr_name2id(_indel_hdr, reg);
    _snp_rid = bcf_hdr_name2id(_snp_hdr, reg);
}

void VcfWriter::setup_vcf_header(
        const char *sn,
        int argc,
        char **argv,
        const char *version,
        const char *ref,
        regions_list_t &regions,
        bcf_hdr_t *hdr,
        bool snp)
{
    time_t today = time(nullptr);
    tm now;
    localtime_r(&today, &now);

    bcf_hdr_printf(hdr, "##fileDate=%d%02d%02d\n", 1900 + now.tm_year, 1 + now.tm_mon, now.tm_mday);
    bcf_hdr_printf(hdr, "##source=%s v%s\n", *argv, version);
    bcf_hdr_printf(hdr, "##reference=%s\n", ref);

    // command
    const size_t BUFFER_SIZE = 0x1000;
    char *buffer = new char[BUFFER_SIZE];

    std::strncpy(buffer, "##command=", BUFFER_SIZE);
    for (int i = 0; i < argc; ++i) {
        std::strncat(buffer, argv[i], BUFFER_SIZE - 1);
        if (i != argc - 1) {
            std::strncat(buffer, " ", 2);
        }
    }
    std::strncat(buffer, "\n", 2);
    bcf_hdr_append(hdr, buffer);

    delete[] buffer;

    // FOMRAT
    bcf_hdr_append(hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
    bcf_hdr_append(hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    bcf_hdr_append(hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled Genotype Likelihood Scores\">\n");
    bcf_hdr_append(hdr, "##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"Reference Read Depth\">\n");
    bcf_hdr_append(hdr, "##FORMAT=<ID=VR,Number=1,Type=Integer,Description=\"Major Variant Read Depth\">\n");

    /*
    if (_opts->report_read_end_score) {
        bcf_hdr_append(hdr, "##FORMAT=<ID=RRES,Number=1,Type=Float,Description=\"Read end score (match)\">\n");
        bcf_hdr_append(hdr, "##FORMAT=<ID=VRES,Number=1,Type=Float,Description=\"Read end score (insertion/deletion)\">\n");
    }
    */

    if (_opts->gvcfoutput) {
        bcf_hdr_append(hdr, "##FORMAT=<ID=DPX,Number=4,Type=Float,Description=\"Minimum, maximum, average, and standard deviation of DP values within this non-variant block\">\n");
        bcf_hdr_append(hdr, "##FORMAT=<ID=RRX,Number=4,Type=Float,Description=\"Minimum, maximum, average, and standard deviation of RR values within this non-variant block\">\n");
        bcf_hdr_append(hdr, "##FORMAT=<ID=VRX,Number=4,Type=Float,Description=\"Minimum, maximum, average, and standard deviation of VR values within this non-variant block\">\n");
    }

    // FILTER
    if (snp) {
        bcf_hdr_append(hdr, "##FILTER=<ID=low_snpqual,Description=\"SNP logistic regression P-value is less than 0.5\">\n");
        bcf_hdr_printf(hdr, "##FILTER=<ID=low_VariantReads,Description=\"Variant read depth is less than %u\">\n", _opts->snp_min_vr);
        //bcf_hdr_printf(hdr, "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than %.2f\">\n", _opts->snp_het_min);
        bcf_hdr_append(hdr, "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than PL cutoff\">\n");
        bcf_hdr_printf(hdr, "##FILTER=<ID=low_coverage,Description=\"Total coverage is less than %u\">\n", _opts->snp_min_dp);
    } else {
        bcf_hdr_append(hdr, "##FILTER=<ID=low_qual,Description=\"Indel logistic regression P-value is less than 0.5\">\n");
        bcf_hdr_printf(hdr, "##FILTER=<ID=low_VariantReads,Description=\"Variant read depth is less than %u\">\n", _opts->indel_min_var_reads);
        bcf_hdr_printf(hdr, "##FILTER=<ID=low_VariantRatio,Description=\"Variant read ratio is less than %.2f\">\n", _opts->indel_min_var_ratio);
        bcf_hdr_printf(hdr, "##FILTER=<ID=low_coverage,Description=\"Total coverage is less than %u\">\n", _opts->indel_min_depth);
    }

    bcf_hdr_printf(hdr, "##FILTER=<ID=high_coverage,Description=\"Total coverage is greater than %u\">\n", _opts->high_cov_cutoff);

    if (snp) {
        bcf_hdr_append(hdr, "##FILTER=<ID=single_strand,Description=\"All variant reads are in a single strand direction\">\n");
    } else {
        bcf_hdr_printf(hdr, "##FILTER=<ID=read_end_ratio,Description=\"Ratio of variants reads within %d bp of read end is greater than %.2f\">\n", _opts->snp_near_end_bases + 2, _opts->indel_max_near_read_end_ratio);
    }

    if (_opts->gvcfoutput) {
        bcf_hdr_append(hdr, "##FILTER=<ID=VRFromDeletion,Description=\"Variant read coverage recorded in non-variant block originates from deletion called upstream\">\n");
    }

    bcf_hdr_append(hdr, "##FILTER=<ID=No_data,Description=\"No valid reads on this site\">\n");
    bcf_hdr_append(hdr, "##FILTER=<ID=No_var,Description=\"No valid variants reads on this site\">\n");

    // INFO
    if (snp) {
        bcf_hdr_append(hdr, "##INFO=<ID=P,Number=1,Type=Float,Description=\"SNP p-value\">\n");
    } else {
        bcf_hdr_append(hdr, "##INFO=<ID=P,Number=1,Type=Float,Description=\"Indel p-value\">\n");
    }

    if (snp) {
        bcf_hdr_append(hdr, "##INFO=<ID=equal_majority,Number=0,Type=Flag,Description=\"The called SNP has an equal number of reads indicating another variant call and base was chosen by highest summed quality score\">\n");
    }

    if (_opts->gvcfoutput) {
        bcf_hdr_printf(hdr, "##INFO=<ID=%s,Number=0,Type=Flag,Description=\"Non-variant block compression scheme %s\">\n", _opts->block_label, _opts->block_label);
        bcf_hdr_append(hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of this non-variant block\">\n");
    }

    // contig
    for (const auto &contig : regions) {
        bcf_hdr_printf(hdr, "##contig=<ID=%s,length=%u>\n", contig.first.c_str(), contig.second);
    }

    bcf_hdr_add_sample(hdr, sn);
    bcf_hdr_add_sample(hdr, nullptr);
}

void VcfWriter::setup_filter_idx()
{
    /**
     * Memoize these so bcf_hdr_id2int doesn't have to be called every single
     * time a variant is called.
     * NOTE this must maintained to match the order and presence of filter
     * names in snp_filter_idx_e and indel_filter_idx_e.
     */
    _snp_filter_idx.push_back(bcf_hdr_id2int(_snp_hdr, BCF_DT_ID, "PASS"));
    _snp_filter_idx.push_back(bcf_hdr_id2int(_snp_hdr, BCF_DT_ID, "No_data"));
    _snp_filter_idx.push_back(bcf_hdr_id2int(_snp_hdr, BCF_DT_ID, "No_var"));
    _snp_filter_idx.push_back(bcf_hdr_id2int(_snp_hdr, BCF_DT_ID, "low_coverage"));
    _snp_filter_idx.push_back(bcf_hdr_id2int(_snp_hdr, BCF_DT_ID, "high_coverage"));
    _snp_filter_idx.push_back(bcf_hdr_id2int(_snp_hdr, BCF_DT_ID, "low_snpqual"));
    _snp_filter_idx.push_back(bcf_hdr_id2int(_snp_hdr, BCF_DT_ID, "low_VariantReads"));
    _snp_filter_idx.push_back(bcf_hdr_id2int(_snp_hdr, BCF_DT_ID, "low_VariantRatio"));
    _snp_filter_idx.push_back(bcf_hdr_id2int(_snp_hdr, BCF_DT_ID, "single_strand"));

    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "PASS"));
    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "No_data"));
    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "No_var"));
    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "low_coverage"));
    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "high_coverage"));
    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "low_qual"));
    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "low_VariantReads"));
    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "low_VariantRatio"));
    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "read_end_ratio"));
    _indel_filter_idx.push_back(bcf_hdr_id2int(_indel_hdr, BCF_DT_ID, "VRFromDeletion"));
}

void VcfWriter::setup_vcf(
        const char *sn,
        int argc,
        char **argv,
        const char *version,
        const char *ref,
        regions_list_t &regions)
{
    this->setup_vcf_header(sn, argc, argv, version, ref, regions, _indel_hdr, false);
    this->setup_vcf_header(sn, argc, argv, version, ref, regions, _snp_hdr, true);
    this->setup_filter_idx();
}

/**
 * Set PL values fro SNPs
 */
void VcfWriter::set_pl(int32_t pl_arr[3], double r, double v)
{
    const double error_rate = 0.001351248680180233;
    const double var_rate = 0.00333333333333;
    const double het_hom_ratio = 2.0;
    const double log_e = std::log10(error_rate);
    const double log_1_e = std::log10(1.0 - error_rate);

    double hom_prob = 1.0 / (1.0 + het_hom_ratio);
    double het_prob = 1.0 - hom_prob;
    double ref_prior = 1.0 - var_rate;
    double hom_prior = var_rate * hom_prob;
    double het_prior = var_rate * het_prob;
    double g[3];
    std::vector< int32_t > g_vec;

    g[0] = (int32_t)std::round(-10.0 * (std::log10(hom_prior) + r * log_e + v * log_1_e));   //  1/1
    g[1] = (int32_t)std::round(-10.0 * (std::log10(het_prior) + (r + v) * std::log10(0.5))); //  0/1
    g[2] = (int32_t)std::round(-10.0 * (std::log10(ref_prior) + r * log_1_e + v * log_e));   //  0/0

    g_vec.push_back(g[0]);
    g_vec.push_back(g[1]);
    g_vec.push_back(g[2]);

    std::sort(g_vec.begin(), g_vec.end());

    pl_arr[2] = g[0] - g_vec[0];
    pl_arr[1] = g[1] - g_vec[0];
    pl_arr[0] = g[2] - g_vec[0];
}

/**
 * Output non-variant gVCF entries for a contiguous stretch of non-variant
 * positions in [last_call_pos, next_var_pos)
 */
void VcfWriter::print_gvcf_span(
        int32_t next_var_pos,
        int32_t last_call_pos,
        int32_t region_end_pos,
        bool is_snp,
        bool prev_del,
        covs_idx_e i_vr,
        covs_idx_e i_rr,
        covs_idx_e i_dp)
{
    if (last_call_pos != _opts->last_call_max) {
        // blegh
        htsFile *fp;
        bcf_hdr_t *hdr;
        bcf1_t *rec;
        int rid, pass_int;
        char *alleles;
        if (is_snp) {
            fp = _snp_fp;
            hdr = _snp_hdr;
            rec = _snp_gvcf_rec;
            rid = _snp_rid;
            pass_int = _snp_pass_int;
            alleles = _snp_alleles_gvcf;
        } else {
            fp = _indel_fp;
            hdr = _indel_hdr;
            rec = _indel_gvcf_rec;
            rid = _indel_rid;
            pass_int = _indel_pass_int;
            alleles = _indel_alleles_gvcf;
        }

        float float_arr[4];
        int32_t tmpi, gt_arr[2], gt_zero, gq;
        int32_t pl_arr[3];
        double block_rr_mean, block_vr_mean;

        gq = bcf_int32_missing;
        gt_zero = bcf_gt_unphased(0);

        GvcfBlock block_vr(_opts->block_abs_lim, _opts->block_rel_lim);
        GvcfBlock block_rr(_opts->block_abs_lim, _opts->block_rel_lim);
        GvcfBlock block_dp(_opts->block_abs_lim, _opts->block_rel_lim);
        int32_t prev_block = last_call_pos;
        bool endofregion = false;
        bool endvariant = false;

        int32_t end_pos = next_var_pos < region_end_pos
            ? next_var_pos
            : region_end_pos;

        coverage_t vr_cov, rr_cov, dp_cov;
        uint32_t k = 0;
        int32_t i = last_call_pos;

        // first pass: k == 0, don't break block and jump to skip_break
        dp_cov = _coverages->get_coverage(i, i_dp);
        rr_cov = _coverages->get_coverage(i, i_rr);
        vr_cov = _coverages->get_coverage(i, i_vr);
        goto skip_first_pass;

        while (i <= end_pos) {
            if (i == next_var_pos) {
                if (i == region_end_pos) {
                    endvariant = true;
                    endofregion = true;
                }
                goto breakblock;
            } else if (i == region_end_pos) {
                endofregion = true;
                goto breakblock;
            }

            vr_cov = _coverages->get_coverage(i, i_vr);
            rr_cov = _coverages->get_coverage(i, i_rr);
            dp_cov = _coverages->get_coverage(i, i_dp);

            if (dp_cov == 0) {
                if (block_dp._min > 0) {
                    goto breakblock;
                }
            } else if (block_dp._max == 0 ||
                       block_dp.break_block((double)dp_cov) ||
                       block_rr.break_block((double)rr_cov) ||
                       block_vr.break_block((double)vr_cov))
            {
                goto breakblock;
            }

            goto skip_break;

        breakblock:
            // break this non-variant block

            bcf_clear1(rec);

            bcf_update_filter(hdr, rec, &pass_int, 1);
            rec->rid = rid;
            rec->pos = prev_block;
            alleles[0] = _sequences->_seq[prev_block];
            bcf_update_alleles_str(hdr, rec, alleles);
            if (block_dp._max == 0) {
                gt_arr[0] = bcf_gt_missing;
                gt_arr[1] = bcf_gt_missing;
            } else if (prev_del && i <= _prev_del_end) {
                gt_arr[0] = _prev_del_gt[0];
                gt_arr[1] = _prev_del_gt[1];
                bcf_add_filter(hdr, rec, _indel_filter_idx[INDEL_FILTER_VRFROMDELETION]);
            } else {
                gt_arr[0] = gt_zero;
                gt_arr[1] = gt_zero;
            }
            bcf_update_genotypes(hdr, rec, gt_arr, 2);

            tmpi = endofregion && !endvariant ? i + 1 : i;
            bcf_update_info_int32(hdr, rec, "END", &tmpi, 1);
            bcf_update_info_flag(hdr, rec, _gvcf_label, nullptr, 1);

            tmpi = (int32_t)block_vr._min;
            bcf_update_format_int32(hdr, rec, "VR", &tmpi, 1);
            tmpi = (int32_t)block_rr._min;
            bcf_update_format_int32(hdr, rec, "RR", &tmpi, 1);
            tmpi = (int32_t)block_dp._min;
            bcf_update_format_int32(hdr, rec, "DP", &tmpi, 1);
            bcf_update_format_int32(hdr, rec, "GQ", &gq, 1);

            block_rr_mean = block_rr.get_mean(k);
            block_vr_mean = block_vr.get_mean(k);
            this->set_pl(pl_arr, (double)block_rr_mean, (double)block_vr_mean);
            bcf_update_format_int32(hdr, rec, "PL", &pl_arr, 3);

            float_arr[0] = (float)block_vr._min;
            float_arr[1] = (float)block_vr._max;
            float_arr[2] = (float)block_vr_mean;
            float_arr[3] = (float)block_vr.get_std_dev(k);
            bcf_update_format_float(hdr, rec, "VRX", float_arr, 4);
            float_arr[0] = (float)block_rr._min;
            float_arr[1] = (float)block_rr._max;
            float_arr[2] = (float)block_rr_mean;
            float_arr[3] = (float)block_rr.get_std_dev(k);
            bcf_update_format_float(hdr, rec, "RRX", float_arr, 4);
            float_arr[0] = (float)block_dp._min;
            float_arr[1] = (float)block_dp._max;
            float_arr[2] = (float)block_dp.get_mean(k);
            float_arr[3] = (float)block_dp.get_std_dev(k);
            bcf_update_format_float(hdr, rec, "DPX", float_arr, 4);

            if (bcf_write1(fp, hdr, rec)) {
                std::cerr << "Failed to write VCF record" << std::endl;
                exit(EXIT_FAILURE);
            }

            if (endofregion) {
                break;
            }

            block_vr.reset();
            block_rr.reset();
            block_dp.reset();
            k = 0;
            prev_block = i;

        skip_break:
            ++i;
        skip_first_pass:
            ++k;
            block_vr.add_value(vr_cov);
            block_rr.add_value(rr_cov);
            block_dp.add_value(dp_cov);
        }
    }
}

void VcfWriter::set_indel_logit_params(
        indel_logit_values_s *values,
        IndelEvent &id_it,
        coverage_t total_coverage)
{
    int32_t var_start = id_it._var_start;
    if (!id_it._isdel && var_start > 0) {
        --var_start;
    }

    double var_ratio = (double)id_it._read_count / (double)total_coverage;

    values->ratio_score = (1.0 - 1.0 / (double)total_coverage) * std::max(var_ratio, 1 - std::abs(std::fma(2, var_ratio, -1)));
    values->strand_dir = id_it.strand_test() ? 1.0 : 0.0;
    values->mean_avnqs = id_it.get_mean_avg_nqs();
    values->seq_entropy = id_it.local_entropy(*_sequences);
    values->mean_var_rate = id_it.get_mean_var_rate();
}

void VcfWriter::set_snp_logit_params(
        snp_logit_values_s *values,
        SnpEvent &sv_it,
        char refbase,
        coverage_t total_coverage)
{
    double var_ratio = (double)sv_it._read_count / (double)total_coverage;

    values->ratio_score = (1.0 - 1.0 / (double)total_coverage) * std::max(var_ratio, 1 - std::abs(std::fma(2, var_ratio, -1)));
    values->base_qual = (double)sv_it.get_qual();
    values->mean_avnqs = sv_it.get_nqs();
    values->rel_pos = sv_it.get_rel_pos();
    values->titv = sv_it.titv(refbase);
}

/**
 * Print SNPs
 */
void VcfWriter::print_snp_buffer(int32_t next_var_pos, const bed_coord_t &seg)
{
    bool equal_majority;
    char refbase;
    double p_value;
    float tmpf;
    int32_t tmpi, gt_arr[2], pl_arr[3];
    coverage_t var_cov, dp_cov, rr_cov, ar_cov, total_coverage;
    qual_t qual;
    snp_logit_values_s logit_values;

    // for each snp event up to next variant or end
    SnpMap::iterator sv_it;
    for (sv_it = _events->_snps.begin(); sv_it != _events->_snps.end() && sv_it->first < next_var_pos; ++sv_it) {
        ar_cov = 0;
        for (auto &allele_it : sv_it->second) {
            size_t num_allele = allele_it.second._read_count;
            ar_cov += num_allele;
        }

        refbase = _sequences->_seq[sv_it->first];
        total_coverage = _coverages->get_coverage(sv_it->first, COVSIDX_SNP_DP) + ar_cov;
        std::memset(pl_arr, 0, sizeof pl_arr);

        // dump snp logit varible values for retraining
        if (_opts->dumpsnp) {
            for (auto &it : sv_it->second) {
                this->set_snp_logit_params(&logit_values, it.second, refbase, total_coverage);

                std::cerr << ">\t" << _region << "," << sv_it->first + 1 << "," << refbase << "," << it.first << "\t"
                          << logit_values.ratio_score << "\t"
                          << logit_values.base_qual << "\t"
                          << logit_values.mean_avnqs << "\t"
                          << logit_values.rel_pos << "\t"
                          << logit_values.titv << std::endl;
            }
        }

        // if variant within a called region
        if (_opts->capturebed == nullptr || (sv_it->first >= seg.first && sv_it->first <= seg.second))
        {
            // coverages
            _snp_counts.clear();
            _snp_prs.clear();

            // sum counts and base quality scores per allele
            for (auto &allele_it : sv_it->second) {
                size_t num_allele = allele_it.second._read_count;
                _snp_counts[allele_it.first] = num_allele;
            }

            // record allele with most supporting reads and its base quality score sum
            char high_base = '\0';
            coverage_t high_base_cov = 0;
            equal_majority = false;

            for (auto &allele_it : sv_it->second) {
                if (_snp_counts[allele_it.first] > high_base_cov) {
                    high_base = allele_it.first;
                    high_base_cov = _snp_counts[high_base];
                    equal_majority = false;
                } else if (_snp_counts[allele_it.first] == high_base_cov) {
                    equal_majority = true;
                }
            }

            // select equal majority snp with highest base quality sum
            if (equal_majority) {
                for (auto &allele_it : sv_it->second) {
                    this->set_snp_logit_params(&logit_values, allele_it.second, refbase, total_coverage);
                    _snp_prs[allele_it.first] = snp_logit(_logit_params->snp, &logit_values);
                }

                double high_pr = _snp_prs[high_base];

                for (auto &allele_it : sv_it->second) {
                    if (_snp_counts[allele_it.first] == high_base_cov && _snp_prs[allele_it.first] > high_pr) {
                        high_pr = _snp_prs[allele_it.first];
                    }
                }
            }

            // P-value and coverage for allele with most supporting reads
            auto &called_snp = sv_it->second.at(high_base);

            var_cov = (coverage_t)std::min(called_snp._read_count, (uint32_t)_coverages->get_max_coverage());
            _coverages->set_coverage(sv_it->first, COVSIDX_SNP_VR, var_cov);
            _coverages->add_coverage(sv_it->first, COVSIDX_SNP_DP, var_cov);

            if (_opts->dumpsnp || _opts->dumpindel) {
                continue;
            }

            this->set_snp_logit_params(&logit_values, called_snp, refbase, total_coverage);
            p_value = snp_logit(_logit_params->snp, &logit_values);

            if (p_value < _opts->min_pr) {
                continue;
            }

            bcf_clear1(_snp_rec);

            rr_cov = _coverages->get_coverage(sv_it->first, COVSIDX_SNP_RR);
            dp_cov = var_cov + rr_cov;
            qual = (qual_t)std::round(-10.0 * std::log10(1.0 - p_value));

            // set filter and genotype
            bcf_update_filter(_snp_hdr, _snp_rec, &_snp_pass_int, 1);
            if (p_value < 0.5) {
                bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_LOW_SNPQUAL]);
            }

            if (dp_cov == 0) {
                bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_NO_DATA]);
                gt_arr[0] = bcf_gt_missing;
                gt_arr[1] = bcf_gt_missing;
            } else {
                if (dp_cov < _opts->snp_min_dp) {
                    bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_LOW_COVERAGE]);
                } else if (dp_cov > _opts->high_cov_cutoff) {
                    bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_HIGH_COVERAGE]);
                }

                if (_opts->enable_snp_strand_cutoff &&
                    dp_cov >= _opts->snp_strand_cutoff &&
                    std::min((double)called_snp._pos_strand / var_cov, (double)(called_snp._read_count - called_snp._pos_strand) / var_cov) < _opts->snp_strand_ratio_cutoff)
                {
                    bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_SINGLE_STRAND]);
                }

                if (var_cov == 0) {
                    bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_NO_VAR]);
                } else if (var_cov <= _opts->snp_min_vr) {
                    bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_LOW_VARIANTREADS]);
                }

                this->set_pl(pl_arr, (double)rr_cov, (double)var_cov);

                if (pl_arr[2] < pl_arr[0]) {
                    if (pl_arr[2] < pl_arr[1]) {
                        gt_arr[0] = bcf_gt_unphased(1);
                        gt_arr[1] = bcf_gt_unphased(1);
                    } else {
                        gt_arr[0] = bcf_gt_unphased(0);
                        gt_arr[1] = bcf_gt_unphased(1);
                    }
                } else {
                    if (pl_arr[1] < pl_arr[0]) {
                        gt_arr[0] = bcf_gt_unphased(0);
                        gt_arr[1] = bcf_gt_unphased(1);
                    } else {
                        gt_arr[0] = bcf_gt_unphased(0);
                        gt_arr[1] = bcf_gt_unphased(0);
                        bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_LOW_VARIANTRATIO]);
                    }
                }
                /*
                double cov_ratio = (double)var_cov / dp_cov;
                if (cov_ratio <= _opts->snp_het_min)
                {
                    gt_arr[0] = bcf_gt_unphased(0);
                    gt_arr[1] = bcf_gt_unphased(0);
                    bcf_add_filter(_snp_hdr, _snp_rec, _snp_filter_idx[SNP_FILTER_LOW_VARIANTRATIO]);
                }
                else if (cov_ratio < _opts->snp_het_max)
                {
                    gt_arr[0] = bcf_gt_unphased(0);
                    gt_arr[1] = bcf_gt_unphased(1);
                } else {
                    gt_arr[0] = bcf_gt_unphased(1);
                    gt_arr[1] = bcf_gt_unphased(1);
                }
                */
            }

            // write gvcf non-variant entries before this variant
            if (_opts->gvcfoutput && _events->_last_call_snp < sv_it->first) {
                this->print_gvcf_span(sv_it->first,
                                    _events->_last_call_snp,
                                    seg.second,
                                    true,
                                    false,
                                    COVSIDX_SNP_VR,
                                    COVSIDX_SNP_RR,
                                    COVSIDX_SNP_DP);
            }
            _events->_last_call_snp = sv_it->first + 1;

            // write variant
            _snp_rec->rid = _snp_rid;
            _snp_rec->pos = sv_it->first;
            _snp_alleles[0] = refbase;
            _snp_alleles[2] = called_snp._allele_base;
            bcf_update_alleles_str(_snp_hdr, _snp_rec, _snp_alleles);
            _snp_rec->qual = qual;
            bcf_update_genotypes(_snp_hdr, _snp_rec, gt_arr, 2);

            tmpf = (float)p_value;
            bcf_update_info_float(_snp_hdr, _snp_rec, "P", &tmpf, 1);

            if (equal_majority) {
                bcf_update_info_flag(_snp_hdr, _snp_rec, "equal_majority", nullptr, 1);
            }

            tmpi = (int32_t)var_cov;
            bcf_update_format_int32(_snp_hdr, _snp_rec, "VR", &tmpi, 1);
            tmpi = (int32_t)rr_cov;
            bcf_update_format_int32(_snp_hdr, _snp_rec, "RR", &tmpi, 1);
            tmpi = (int32_t)dp_cov;
            bcf_update_format_int32(_snp_hdr, _snp_rec, "DP", &tmpi, 1);
            bcf_update_format_int32(_snp_hdr, _snp_rec, "GQ", &_gt_missing, 1);
            bcf_update_format_int32(_snp_hdr, _snp_rec, "PL", &pl_arr, 3);

            if (bcf_write1(_snp_fp, _snp_hdr, _snp_rec)) {
                std::cerr << "Failed to write VCF record" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    // Clear processed SNPs from _snps
    _events->_snps.erase(_events->_snps.begin(), sv_it);
}

/**
 * Process indels
 */
void VcfWriter::print_indel_buffer(int32_t next_var_pos, const bed_coord_t &seg)
{
    char refbase;
    double p_value, var_ratio, ref_ratio;
    float tmpf;
    int32_t tmpi, var_start, gt_arr[2], pl_arr[3];
    std::string alt, ref;
    coverage_t ar_cov, total_coverage, total_ref_cov;
    qual_t qual;
    indel_logit_values_s logit_values;

    IndelMap::iterator iv_it;
    for (iv_it = _events->_indels.begin(); iv_it != _events->_indels.end() && iv_it->first < next_var_pos; ++iv_it) {
        // dump indel logit varible values for retraining
        if (_opts->dumpindel) {
            for (auto id_it : iv_it->second) {
                var_start = id_it.second._var_start;
                if (!id_it.second._isdel && var_start > 0) {
                    --var_start;
                }

                if ((total_coverage = _coverages->get_coverage(var_start, COVSIDX_INDEL_DP)) == 0) {
                    // Case where only supporting reads are unanchored insertions
                    continue;
                }

                this->set_indel_logit_params(&logit_values, id_it.second, total_coverage);

                ref.clear();
                alt.clear();

                if (id_it.second._isdel) {
                    ref.append(_sequences->_seq + id_it.second._var_start - 1, (size_t)id_it.second._var_len + 1);
                    alt.push_back(_sequences->_seq[id_it.second._var_start - 1]);
                } else {
                    refbase = _sequences->_seq[id_it.second._var_start - 1];
                    ref.push_back(refbase);
                    alt.push_back(refbase);
                    alt.append(id_it.second._seq);
                }

                std::cerr << ">\t" << _region << "," << iv_it->first << "," << ref << "," << alt << "\t"
                          << logit_values.ratio_score << "\t"
                          << logit_values.strand_dir << "\t"
                          << logit_values.mean_avnqs << "\t"
                          << logit_values.seq_entropy << "\t"
                          << logit_values.mean_var_rate << std::endl;
            }
        }

        ar_cov = 0;
        auto ivg_it = iv_it->second.begin();
        IndelEvent *called_indel_ptr = &ivg_it->second;

        auto called_indel_read_count = (coverage_t)std::min(called_indel_ptr->_read_count, (uint32_t)_coverages->get_max_coverage());

        // find variant with maximum supporting reads
        while (ivg_it != iv_it->second.end()) {
            ar_cov += ivg_it->second._read_count;
            if (ivg_it->second._read_count > called_indel_read_count) {
                called_indel_ptr = &ivg_it->second;
            }
            ++ivg_it;
        }

        IndelEvent &called_indel = *called_indel_ptr;

        // if variant within a called region
        if (_opts->capturebed == nullptr ||
            (called_indel._var_start >= seg.first &&
             called_indel._var_start <= seg.second))
        {
            // "start position" for coverage
            var_start = called_indel._var_start;
            // what happens if var_start == 0, is that even valid?
            if (!called_indel._isdel && var_start > 0) {
                --var_start;
            }

            if ((total_coverage = _coverages->get_coverage(var_start, COVSIDX_INDEL_DP)) == 0) {
                // Case where only supporting reads are unanchored insertions
                continue;
            }

            this->set_indel_logit_params(&logit_values, called_indel, total_coverage);
            p_value = indel_logit(_logit_params->indel, &logit_values);
            qual = (qual_t)std::round(-10.0 * std::log10(1.0 - p_value));

            // record coverage for deletion span
            if (called_indel._isdel) {
                _coverages->add_coverage_range(var_start, COVSIDX_INDEL_VR, called_indel._var_len, called_indel_read_count);
            }

            if (_opts->dumpsnp || _opts->dumpindel || p_value < _opts->min_pr) {
                continue;
            }

            bcf_clear1(_indel_rec);
            bcf_update_filter(_indel_hdr, _indel_rec, &_indel_pass_int, 1);

            // set filter
            var_ratio = (double)called_indel_read_count / total_coverage;

            if (called_indel_read_count < _opts->indel_min_var_reads) {
                bcf_add_filter(_indel_hdr, _indel_rec, _indel_filter_idx[INDEL_FILTER_LOW_VARIANTREADS]);
            }

            if (total_coverage < _opts->indel_min_depth) {
                bcf_add_filter(_indel_hdr, _indel_rec, _indel_filter_idx[INDEL_FILTER_LOW_COVERAGE]);
            } else if (total_coverage > _opts->high_cov_cutoff) {
                bcf_add_filter(_indel_hdr, _indel_rec, _indel_filter_idx[INDEL_FILTER_HIGH_COVERAGE]);
            }

            if (var_ratio < _opts->indel_min_var_ratio) {
                bcf_add_filter(_indel_hdr, _indel_rec, _indel_filter_idx[INDEL_FILTER_LOW_VARIANTRATIO]);
            }
            /*
            if (_opts->indel_strand_dir_filter && !called_indel.strand_test()) {
                bcf_add_filter(_indel_hdr, _indel_rec, _indel_filter_idx[INDEL_FILTER_SINGLE_STRAND]);
            }
            */
            if (called_indel.get_near_read_end_ratio() > _opts->indel_max_near_read_end_ratio) {
                bcf_add_filter(_indel_hdr, _indel_rec, _indel_filter_idx[INDEL_FILTER_READ_END_RATIO]);
            }

            // total coverage
            total_ref_cov = _coverages->get_coverage(var_start, COVSIDX_INDEL_RR);
            if (!called_indel._isdel)
                total_ref_cov -= _coverages->get_coverage(var_start, COVSIDX_INDEL_RR_INS);

            this->set_pl(pl_arr, (double)total_ref_cov, (double)called_indel_read_count);

            /* Use hard cutoffs for now until we get better priors for PL values
            if (pl_arr[2] < pl_arr[0]) {
                if (pl_arr[2] < pl_arr[1]) {
                    gt_arr[0] = bcf_gt_unphased(1);
                    gt_arr[1] = bcf_gt_unphased(1);
                }
                else {
                    gt_arr[0] = bcf_gt_unphased(0);
                    gt_arr[1] = bcf_gt_unphased(1);
                }
            } else {
                if (pl_arr[1] < pl_arr[0]) {
                    gt_arr[0] = bcf_gt_unphased(0);
                    gt_arr[1] = bcf_gt_unphased(1);
                } else {
                    gt_arr[0] = bcf_gt_unphased(0);
                    gt_arr[1] = bcf_gt_unphased(0);
                }
            }
            */

            double score, read_end_score[3];
            double near_read_end_sum[3];
            uint32_t near_read_end_n[3];
            for (size_t i = 0; i < 3; ++i) {
                near_read_end_sum[i] = 0.0;
                near_read_end_n[i] = 0;
            }

            // near read end score
            // for each read
            for (auto &bc_it : *_cigar_list) {
                if (bc_it.first.first > called_indel._var_start) {
                    break;
                }
                if (bc_it.first.second < called_indel._var_start) {
                    continue;
                }

                int32_t r_end = bc_it.first.first;
                int32_t q_end = bc_it.first.first;

                // for each cigar
                for (auto &cigar_it : bc_it.second) {
                    if (bam_cigar_type(cigar_it) & 0x01) {
                        q_end += bam_cigar_oplen(cigar_it);
                    }
                    if (bam_cigar_type(cigar_it) & 0x02) {
                        r_end += bam_cigar_oplen(cigar_it);
                    }

                    if (q_end > called_indel._var_start || r_end > called_indel._var_start) {
                        size_t idx;
                        switch (bam_cigar_op(cigar_it)) {
                        case BAM_CMATCH: idx = 0; break;
                        case BAM_CINS:   idx = 1; break;
                        case BAM_CDEL:   idx = 2; break;
                        default:         idx = 3; break;
                        }

                        if (idx > 2) {
                            continue;
                        }

                        if (bam_cigar_op(cigar_it) == BAM_CMATCH) /*||
                            (_opts->report_read_end_score &&
                             (bam_cigar_op(cigar_it) == BAM_CINS ||
                              bam_cigar_op(cigar_it) == BAM_CDEL)))*/
                        {
                            score = (double)(called_indel._var_start - bc_it.first.first) / (double)(bc_it.first.second - bc_it.first.first);
                            if (score > 0.5) {
                                score = 1.0 - score;
                            }
                            near_read_end_sum[idx] += score;
                            ++near_read_end_n[idx];
                        }

                        break;
                    }
                }
            }
            for (size_t i = 0; i < 3; ++i) {
                read_end_score[i] = near_read_end_n[i] != 0
                    ? near_read_end_sum[i] / (double)near_read_end_n[i]
                    : 0.0;
            }

            double n_tmp = (called_indel_read_count == 0) ? total_coverage : called_indel_read_count + total_ref_cov;
            var_ratio = (double)called_indel_read_count / n_tmp;
            ref_ratio = (double)total_ref_cov / n_tmp;
            if (var_ratio <= _opts->indel_het_min) {
                gt_arr[0] = bcf_gt_unphased(0);
                gt_arr[1] = bcf_gt_unphased(0);
            } else if (var_ratio < _opts->indel_het_max) {
                if (ref_ratio * read_end_score[0] < 0.06) {
                    gt_arr[0] = bcf_gt_unphased(1);
                    gt_arr[1] = bcf_gt_unphased(1);
                } else {
                    gt_arr[0] = bcf_gt_unphased(0);
                    gt_arr[1] = bcf_gt_unphased(1);
                }
            } else {
                gt_arr[0] = bcf_gt_unphased(1);
                gt_arr[1] = bcf_gt_unphased(1);
            }

            // ref and alt sequences
            ref.clear();
            alt.clear();

            if (called_indel._isdel) {
                ref.append(_sequences->_seq + called_indel._var_start - 1, (size_t)called_indel._var_len + 1);
                alt.push_back(_sequences->_seq[called_indel._var_start - 1]);
            } else {
                refbase = _sequences->_seq[called_indel._var_start - 1];
                ref.push_back(refbase);
                alt.push_back(refbase);
                alt.append(called_indel._seq);
            }

            if (p_value < 0.5) {
                bcf_add_filter(_indel_hdr, _indel_rec, _indel_filter_idx[INDEL_FILTER_LOW_QUAL]);
            }

            // write gvcf non-variant entries before this variant
            if (_opts->gvcfoutput && _events->_last_call_indel < called_indel._var_start - 1) {
                this->print_gvcf_span(
                    called_indel._var_start - 1,
                    _events->_last_call_indel,
                    seg.second,
                    false,
                    _prev_del,
                    COVSIDX_INDEL_VR,
                    COVSIDX_INDEL_RR,
                    COVSIDX_INDEL_DP
                );
            }
            _events->_last_call_indel = called_indel._var_start;
            if ((_prev_del = called_indel._isdel)) {
                _prev_del_end = var_start + called_indel._var_len;
                _prev_del_gt[0] = gt_arr[0];
                _prev_del_gt[1] = gt_arr[1];
            }

            _indel_rec->rid = _indel_rid;
            _indel_rec->pos = called_indel._var_start - 1;
            _indel_alleles[0] = ref.c_str();
            _indel_alleles[1] = alt.c_str();
            bcf_update_alleles(_indel_hdr, _indel_rec, _indel_alleles, 2);

            _indel_rec->qual = qual;
            bcf_update_genotypes(_indel_hdr, _indel_rec, gt_arr, 2);

            tmpf = (float)p_value;
            bcf_update_info_float(_indel_hdr, _indel_rec, "P", &tmpf, 1);

            tmpi = (int32_t)called_indel_read_count;
            bcf_update_format_int32(_indel_hdr, _indel_rec, "VR", &tmpi, 1);
            tmpi = (int32_t)total_ref_cov;
            bcf_update_format_int32(_indel_hdr, _indel_rec, "RR", &tmpi, 1);
            tmpi = (int32_t)total_coverage;
            bcf_update_format_int32(_indel_hdr, _indel_rec, "DP", &tmpi, 1);
            bcf_update_format_int32(_indel_hdr, _indel_rec, "GQ", &_gt_missing, 1);
            bcf_update_format_int32(_indel_hdr, _indel_rec, "PL", &pl_arr, 3);

            /*
            if (_opts->report_read_end_score) {
                tmpf = (float)read_end_score[0];
                bcf_update_format_float(_indel_hdr, _indel_rec, "RRES", &tmpf, 1);
                tmpf = (float)read_end_score[called_indel._isdel ? 2: 1];
                bcf_update_format_float(_indel_hdr, _indel_rec, "VRES", &tmpf, 1);
            }
            */

            if (bcf_write1(_indel_fp, _indel_hdr, _indel_rec)) {
                std::cerr << "Failed to write VCF record" << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    // Clear processed indels from _indels
    _events->_indels.erase(_events->_indels.begin(), iv_it);
}
