#ifndef _XATLAS_SNPEVENT_H
#define _XATLAS_SNPEVENT_H

#include "xatlas/Xatlas.h"

class SnpEvent
{
  private:
    SnpEvent();

  public:
    char _allele_base;
    uint32_t _read_count;
    uint32_t _pos_strand;
    qual_t _qual;
    double _rel_pos;
    double _nqs;

    SnpEvent(char alt, bool rs, qual_t qual, double rp, double nq);
    SnpEvent(const SnpEvent &);
    ~SnpEvent();
    void add_snp_event(const SnpEvent &sv);
    double get_qual();
    double get_rel_pos();
    double get_nqs();
    double titv(char ref);
};

#endif /* _XATLAS_SNPEVENT_H */
