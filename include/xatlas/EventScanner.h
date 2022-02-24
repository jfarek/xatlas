#ifndef _XATLAS_EVENTHOLDER_H
#define _XATLAS_EVENTHOLDER_H

#include "htslib/sam.h"

#include "xatlas/IndelEvent.h"
#include "xatlas/SnpEvent.h"
#include "xatlas/Xatlas.h"

namespace xatlas {

/* position -> allele/variant id -> variant set */
typedef std::map< int32_t, std::map< char, SnpEvent > > SnpMap;
typedef std::map< int32_t, std::map< std::string, IndelEvent > > IndelMap;
typedef std::vector< std::pair< int32_t, SnpEvent > > AReadsSnpList;
typedef std::vector< IndelEvent > AReadsIndelList;

class EventScanner
{
  private:
    int32_t _near_end;
    double _snp_max_sub;
    double _snp_max_gap;

    EventScanner();
    double snp_nqs(int32_t p, int32_t qlen, const uint8_t *qquals);
    int32_t left_align(const ReferenceSequence &refseq, std::string& seq, int32_t indel_pos, int32_t read_start_pos);

  public:
    SnpMap _snps;
    IndelMap _indels;
    int32_t _last_call_indel;
    int32_t _last_call_snp;

    EventScanner(int32_t near_end, double snp_max_sub, double snp_max_gap);
    EventScanner(const EventScanner &);
    ~EventScanner();
    void reset();
    void collect_indels(
        bam1_t *read,
        const ReferenceSequence &sequences,
        CoverageCounter &coverages,
        bed_coord_t &segment);
    void collect_snps(
        bam1_t *read,
        const ReferenceSequence &sequences,
        CoverageCounter &coverages,
        bed_coord_t &segment);
};

} // namepace xatlas

#endif /* _XATLAS_EVENTHOLDER_H */
