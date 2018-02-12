#include "ReferenceSequence.hpp"
#include <cstring>

/**
 * Handle reference FASTA
 */

ReferenceSequence::ReferenceSequence()
    : _seq(nullptr),
      _len(0),
      _status(FAI_OKAY)
{
}

ReferenceSequence::ReferenceSequence(const char *refseq)
{
    if ((_fai = fai_load(refseq)) == nullptr) {
        // failed to load reference index
        _status = FAI_FAILURE;
    } else {
        // loaded reference index successfully
        _status = FAI_OKAY;
        _seq = nullptr;
        _len = 0;
    }
}

ReferenceSequence::ReferenceSequence(const ReferenceSequence &) = default;

ReferenceSequence::~ReferenceSequence()
{
    fai_destroy(_fai);
    std::free(_seq);
}

void ReferenceSequence::set_region(const char *region)
{
    int32_t region_len;

    // free previous region
    std::free(_seq);

    // set current region
    _seq = fai_fetch(_fai, region, &region_len);
    _len = region_len;

    // set all reference sequence to upper case
    for (size_t i = 0; i < (size_t)region_len; ++i) {
        _seq[i] = std::toupper(_seq[i]);
    }
}
