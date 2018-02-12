#ifndef _XATLAS_BAM_H
#define _XATLAS_BAM_H

#include "htslib/faidx.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

enum bam_status {
    BAM_OKAY = 0,
    BAM_BAD_FILE,
    BAM_BAD_HEADER,
    BAM_BAD_INDEX,
    BAM_BAD_REF_INDEX
};
typedef enum bam_status bam_status_e;

class Bam
{
  private:
    hts_idx_t *_idx;
    htsThreadPool *_thread_pool;

    Bam();

  public:
    samFile *_sf;
    hts_itr_t *_iter;
    bam_hdr_t *_hdr;
    bam_status_e _status;

    explicit Bam(const char *sf_fn, const char *ref_fn, uint8_t nthreads);
    Bam(const Bam &);
    ~Bam();
    void set_iter(const char *where);
};

#endif /* _XATLAS_BAM_H */
