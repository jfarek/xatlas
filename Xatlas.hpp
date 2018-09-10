#ifndef _XATLAS_H
#define _XATLAS_H

#include <cstdint>
#include <deque>
#include <map>
#include <string>
#include <vector>

#define XATLAS_NAME    "xAtlas"
#define XATLAS_VERSION "0.2.1"

typedef uint16_t qual_t;
typedef uint16_t coverage_t;
typedef std::pair< int32_t, int32_t > bed_coord_t;
typedef std::vector< bed_coord_t > bed_coord_list_t;
typedef std::map< std::string, bed_coord_list_t > bed_coord_map_t;
typedef std::pair< std::string, bed_coord_list_t > bed_chr_t;
typedef std::vector< bed_chr_t > bed_list_t;
typedef std::vector< std::pair< std::string, uint32_t > > regions_list_t;
typedef std::deque< std::pair< bed_coord_t, std::vector< uint32_t > > > cigar_list_t;

typedef struct opts {
    bool gvcfoutput;
    bool bgzf;
    bool dumpsnp;
    bool dumpindel;
    /*bool report_read_end_score;*/
    const char *capturebed;
    double block_abs_lim;
    double block_rel_lim;
    double block_rel_min;
    char block_label[24];
    int32_t last_call_max;
    coverage_t high_cov_cutoff;
    double min_pr;

    coverage_t snp_min_dp;
    coverage_t snp_min_vr;
    bool enable_snp_strand_cutoff;
    coverage_t snp_strand_cutoff;
    double snp_strand_ratio_cutoff;
    int32_t snp_near_end_bases;
    /*double snp_het_min;*/
    /*double snp_het_max;*/

    coverage_t indel_min_depth;
    coverage_t indel_min_var_reads;
    double indel_min_var_ratio;
    double indel_max_near_read_end_ratio;
    double indel_het_min;
    double indel_het_max;
    /*bool indel_strand_dir_filter;*/

} opts_s;

#endif /* _XATLAS_H */
