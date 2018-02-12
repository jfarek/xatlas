#ifndef _XATLAS_INDELEVENT_H
#define _XATLAS_INDELEVENT_H

#include "CoverageCounter.hpp"
#include "ReferenceSequence.hpp"
#include "Xatlas.hpp"
#include <cstdint>

enum strand_mask {
    strand_mask_none = 0x00,
    strand_mask_pos = 0x01,
    strand_mask_neg = 0x02,
    strand_mask_both = 0x03
};
typedef enum strand_mask strand_mask_e;

class IndelEvent
{
  private:
    uint8_t _strand_mask;
    qual_t _map_qual;

    IndelEvent();

  public:
    int32_t _var_start;
    int32_t _var_len;
    int32_t _q_pos;
    int32_t _r_pos;
    std::string _seq;
    std::string _id;
    uint16_t _read_count;
    bool _isdel;
    uint16_t _near_read_end_count;
    qual_t _avg_nbq;
    double _var_rate_gap_and_mismatch;

    IndelEvent(
        int32_t start,
        int32_t len,
        int32_t q_pos,
        int32_t r_pos,
        std::string &seq,
        bool strand,
        qual_t mqual,
        uint16_t nearend = 0,
        qual_t avgnbq = 0,
        double vrate = 0.0);
    IndelEvent(const IndelEvent &);
    ~IndelEvent();
    void add_indel_event(const IndelEvent &iv);
    bool strand_test();
    double local_entropy(ReferenceSequence &sequences);
    double get_near_read_end_ratio();
    double get_mean_avg_nqs();
    double get_mean_var_rate();
};

#endif /* _XATLAS_INDELEVENT_H */
