#ifndef _XATLAS_SEQUENCE_H
#define _XATLAS_SEQUENCE_H

#include "Xatlas.hpp"
#include "htslib/faidx.h"

enum fai_status {
    FAI_OKAY = 0,
    FAI_FAILURE
};
typedef enum fai_status fai_status_e;

class ReferenceSequence
{
  private:
    ReferenceSequence();

  public:
    faidx_t *_fai;
    char *_seq;
    int32_t _len;
    fai_status_e _status;

    explicit ReferenceSequence(const char *refseq);
    ReferenceSequence(const ReferenceSequence &);
    ~ReferenceSequence();
    void set_region(const char *region);
};

#endif /* _XATLAS_SEQUENCE_H */
