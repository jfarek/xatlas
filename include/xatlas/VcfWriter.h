#ifndef _XATLAS_VCFWRITER_H
#define _XATLAS_VCFWRITER_H

#include "xatlas/CoverageCounter.h"
#include "xatlas/EventScanner.h"
#include "xatlas/Logit.h"
#include "xatlas/Xatlas.h"
#include "htslib/vcf.h"
#include <cstdint>
#include <string>

enum snp_filter_idx {
    SNP_FILTER_PASS = 0,
    SNP_FILTER_NO_DATA,
    SNP_FILTER_NO_VAR,
    SNP_FILTER_LOW_COVERAGE,
    SNP_FILTER_HIGH_COVERAGE,
    SNP_FILTER_LOW_SNPQUAL,
    SNP_FILTER_LOW_VARIANTREADS,
    SNP_FILTER_LOW_VARIANTRATIO,
    SNP_FILTER_SINGLE_STRAND
};
typedef enum snp_filter_idx snp_filter_idx_e;

enum indel_filter_idx {
    INDEL_FILTER_PASS = 0,
    INDEL_FILTER_NO_DATA,
    INDEL_FILTER_NO_VAR,
    INDEL_FILTER_LOW_COVERAGE,
    INDEL_FILTER_HIGH_COVERAGE,
    INDEL_FILTER_LOW_QUAL,
    INDEL_FILTER_LOW_VARIANTREADS,
    INDEL_FILTER_LOW_VARIANTRATIO,
    INDEL_FILTER_READ_END_RATIO,
    INDEL_FILTER_VRFROMDELETION
    /*INDEL_FILTER_SINGLE_STRAND*/
};
typedef enum snp_filter_idx snp_filter_idx_e;

class VcfWriter
{
  private:
    CoverageCounter *_coverages;
    ReferenceSequence *_sequences;
    EventScanner *_events;
    cigar_list_t *_cigar_list;
    std::string _region;
    int _indel_rid;
    int _snp_rid;
    int _indel_pass_int;
    int _snp_pass_int;
    bcf_hdr_t *_indel_hdr;
    bcf_hdr_t *_snp_hdr;
    htsFile *_indel_fp;
    htsFile *_snp_fp;
    const opts_s *_opts;
    const logit_params *_logit_params;
    bcf1_t *_indel_rec;
    bcf1_t *_indel_gvcf_rec;
    bcf1_t *_snp_rec;
    bcf1_t *_snp_gvcf_rec;
    const char *_gvcf_label;
    int32_t _prev_del_end;
    int32_t _prev_del_gt[2];
    int _gt_missing;

    char _snp_alleles[4];
    const char *_indel_alleles[2];
    char _snp_alleles_gvcf[4];
    char _indel_alleles_gvcf[4];

    std::map< char, uint32_t > _snp_counts;
    std::map< char, double > _snp_prs;

    std::vector< int > _snp_filter_idx;
    std::vector< int > _indel_filter_idx;

    VcfWriter();
    void setup_filter_idx();
    void setup_vcf_header(
        const char *sn,
        int argc,
        char **argv,
        const char *version,
        const char *ref,
        regions_list_t &regions,
        bcf_hdr_t *hdr,
        bool snp);

  public:
    bool _prev_del;

    VcfWriter(CoverageCounter &coverages,
              ReferenceSequence &sequences,
              EventScanner &events,
              cigar_list_t &cigar_list,
              bcf_hdr_t *indel_hdr,
              bcf_hdr_t *snp_hdr,
              htsFile *indel_fp,
              htsFile *snp_fp,
              const opts_s *opts,
              const logit_params *logit_params);
    VcfWriter(const VcfWriter &);
    ~VcfWriter();
    void set_region(const char *reg);
    void setup_vcf(
        const char *sn,
        int argc,
        char **argv,
        const char *version,
        const char *ref,
        regions_list_t &regions);    void set_pl(int32_t pl_arr[3], double r, double v);
    void print_gvcf_span(
        int32_t next_var_pos,
        int32_t last_call_pos,
        int32_t region_end_pos,
        bool is_snp,
        bool prev_del,
        covs_idx_e i_vr,
        covs_idx_e i_rr,
        covs_idx_e i_dp);
    void set_indel_logit_params(
        indel_logit_values_s *values,
        IndelEvent &id_it,
        coverage_t total_coverage);
    void set_snp_logit_params(
        snp_logit_values_s *values,
        SnpEvent &sv_it,
        char refbase,
        coverage_t total_coverage);
    void print_snp_buffer(int32_t next_var_pos, const bed_coord_t &seg);
    void print_indel_buffer(int32_t next_var_pos, const bed_coord_t &seg);
};

#endif /* _XATLAS_VCFWRITER_H */
