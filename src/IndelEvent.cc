#include <cmath>
#include <cstring>

#include "xatlas/IndelEvent.h"

namespace xatlas {

IndelEvent::IndelEvent()
    : _strand_mask(strand_mask_none),
      _map_qual(0),
      _var_start(0),
      _var_len(0),
      _q_pos(0),
      _seq(""),
      _id(""),
      _read_count(0),
      _isdel(false),
      _near_read_end_count(0),
      _avg_nbq(0),
      _var_rate_gap_and_mismatch(0)
{
}

IndelEvent::IndelEvent(int32_t start,
                       int32_t len,
                       int32_t q_pos,
                       std::string &seq,
                       bool strand,
                       qual_t mqual,
                       uint16_t nearend,
                       qual_t avgnbq,
                       double vrate)
    : _strand_mask(strand ? strand_mask_pos : strand_mask_neg),
      _map_qual(mqual),
      _var_start(start),
      _var_len(len),
      _q_pos(q_pos),
      _seq(seq),
      _id(seq.empty() ? std::to_string(len) : seq),
      _read_count(1),
      _isdel(seq.empty()),
      _near_read_end_count(nearend),
      _avg_nbq(avgnbq),
      _var_rate_gap_and_mismatch(vrate)
{
}

IndelEvent::IndelEvent(const IndelEvent &) = default;

IndelEvent::~IndelEvent() = default;

void IndelEvent::add_indel_event(const IndelEvent &iv)
{
    ++_read_count;
    _strand_mask |= iv._strand_mask;
    _near_read_end_count += iv._near_read_end_count;
    _map_qual += iv._map_qual;
    _avg_nbq += iv._avg_nbq;
    _var_rate_gap_and_mismatch += iv._var_rate_gap_and_mismatch;
}

bool IndelEvent::strand_test()
{
    return !(bool)(~_strand_mask & strand_mask_both);
}

/* TODO figure out wtf this is doing */
double IndelEvent::local_entropy(ReferenceSequence &refseq)
{
    int32_t lower, upper, len;

    if (_isdel) {
        lower = _var_start <= 10 ? 1 : _var_start - 9;
        upper = _var_start + _var_len + 11;
    } else {
        lower = _var_start <= 10 ? 1 : _var_start - 10;
        upper = _var_start + 12;
    }
    len = upper - lower;

    if (len < _var_len) {
        return 0.0;
    }

    bool found, match;
    char *p, *q, *r = refseq._seq + lower - 1;
    uint16_t tmp[44];
    int32_t tmp_pos, j, k, i = 0, n = 0, end = len - _var_len;

    std::memset(tmp, 0, 44 * sizeof(uint16_t));

    while (i < end) {
        found = false;
        p = r + i;
        j = 0;

        while (j < n) {
            match = true;
            tmp_pos = 2 * j;
            q = r + tmp[tmp_pos];
            k = 0;

            while (k < _var_len) {
                if (p[k] != q[k]) {
                    match = false;
                    break;
                }

                ++k;
            }

            if (match) {
                found = true;
                ++tmp[tmp_pos + 1];
                break;
            }

            ++j;
        }

        if (!found) {
            tmp_pos = 2 * n;
            tmp[tmp_pos] = i;
            tmp[tmp_pos + 1] = 1;
            ++n;
        }

        ++i;
    }

    double e, f, s;
    e = 0.0;
    s = 1.0 / (1.0 + end);
    tmp_pos = 2 * n + 1;

    for (j = 1; j < tmp_pos; j += 2) {
        f = s * (double)tmp[j];
        e -= f * std::log(f);
    }

    return e;
}

double IndelEvent::get_near_read_end_ratio()
{
    return (double)_near_read_end_count / _read_count;
}

double IndelEvent::get_mean_avg_nqs()
{
    return (double)_avg_nbq / _read_count;
}

double IndelEvent::get_mean_var_rate()
{
    return (double)_var_rate_gap_and_mismatch / _read_count;
}

} // namepace xatlas
