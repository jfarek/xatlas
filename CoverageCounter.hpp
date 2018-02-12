#ifndef _XATLAS_COVERAGECOUNTER_H
#define _XATLAS_COVERAGECOUNTER_H

#include "Xatlas.hpp"

#define COVS_ALIGN 7

enum covs_idx {
    COVSIDX_INDEL_VR = 0,
    COVSIDX_INDEL_RR = 1,
    COVSIDX_INDEL_RR_INS = 2,
    /*COVSIDX_INDEL_AR,*/
    COVSIDX_INDEL_DP = 3,
    COVSIDX_SNP_VR = 4,
    COVSIDX_SNP_RR = 5,
    /*COVSIDX_SNP_AR,*/
    COVSIDX_SNP_DP = 6
};
typedef covs_idx covs_idx_e;

class CoverageCounter
{
  private:
    uint16_t (*_covs)[COVS_ALIGN];
    uint32_t _len;
    coverage_t _max_cov;

    CoverageCounter();

  public:
    CoverageCounter(size_t len, coverage_t max_cov);
    CoverageCounter(const CoverageCounter &);
    ~CoverageCounter();
    void reset();
    void add_coverage(int32_t pos, covs_idx_e idx);
    void add_coverage(int32_t pos, covs_idx_e idx, coverage_t add_cov);
    void add_coverage_range(int32_t pos, covs_idx_e idx, int32_t len);
    void add_coverage_range(int32_t pos, covs_idx_e idx, int32_t len, coverage_t add_cov);
    void set_coverage(int32_t pos, covs_idx_e idx, coverage_t cov);
    coverage_t get_coverage(int32_t pos, covs_idx_e idx);
    bool test_high_coverage(int32_t pos, covs_idx_e idx, int32_t l_qseq);
};

#endif /* _XATLAS_COVERAGECOUNTER_H */
