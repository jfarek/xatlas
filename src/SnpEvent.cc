#include "xatlas/SnpEvent.h"
#include <cmath>

SnpEvent::SnpEvent()
    : _allele_base(0),
      _read_count(0),
      _pos_strand(0),
      _qual(0),
      _rel_pos(0.0),
      _nqs(0.0)
{
}

SnpEvent::SnpEvent(char alt,
                   bool rs,
                   qual_t qual,
                   double rp,
                   double nq)
    : _allele_base(alt),
      _read_count(0),
      _pos_strand(0),
      _qual(qual),
      _rel_pos(rp),
      _nqs(nq)
{
    ++_read_count;
    if (rs) {
        ++_pos_strand;
    }
}

SnpEvent::SnpEvent(const SnpEvent &) = default;

SnpEvent::~SnpEvent() = default;

void SnpEvent::add_snp_event(const SnpEvent &sv)
{
    _read_count += sv._read_count;
    _pos_strand += sv._pos_strand;
    _qual += sv._qual;
    _rel_pos += sv._rel_pos;
    _nqs += sv._nqs;
}

double SnpEvent::get_qual()
{
    return (double)_qual / (double)_read_count;
}

double SnpEvent::get_rel_pos()
{
    return 1.0 - std::abs(0.5 - (double)_rel_pos / (double)_read_count);
}

double SnpEvent::get_nqs()
{
    return (double)_nqs / (double)_read_count;
}

double SnpEvent::titv(char ref)
{
    // transversion
    double titv = 0.0;

    // transition
    switch (ref) {
    case 'A': if (_allele_base == 'G') { titv = 1.0; } break;
    case 'C': if (_allele_base == 'T') { titv = 1.0; } break;
    case 'G': if (_allele_base == 'A') { titv = 1.0; } break;
    case 'T': if (_allele_base == 'C') { titv = 1.0; } break;
    default: break;
    }

    return titv;
}
