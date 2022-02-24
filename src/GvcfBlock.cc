#include <algorithm>
#include <cmath>

#include "xatlas/GvcfBlock.h"

namespace xatlas {

/**
 * Calculations for non-variant gVCF blocks
 */

GvcfBlock::GvcfBlock()
    : _abs_lim(0.0),
      _rel_lim(0.0),
      _bound(0.0)
{
    this->reset();
}

GvcfBlock::GvcfBlock(double block_abs_lim, double block_rel_lim)
    : _abs_lim(block_abs_lim),
      _rel_lim(block_rel_lim),
      _bound(block_abs_lim / block_rel_lim)
{
    this->reset();
}

GvcfBlock::GvcfBlock(const GvcfBlock &) = default;

GvcfBlock::~GvcfBlock() = default;

/* get*(): each assumes k is always > 0 */

// Arithmetic mean
double GvcfBlock::get_mean(uint32_t k)
{
    return _sum / k;
}

// Sample standard deviation
double GvcfBlock::get_std_dev(uint32_t k)
{
    return (k > 1 && _min != _max)
               ? std::sqrt(std::fma(-_sum, _sum / k, _sum_sq) / (k - 1.0))
               : 0.0;
}

// Would value invalidate current block if added?
bool GvcfBlock::break_block(double value)
{
    return (value > _high ||
            value < _max - (value < _bound ? _abs_lim : _rel_lim * value));
}

// Value must be valid to extend block
void GvcfBlock::add_value(double value)
{
    if (value < _min) {
        _min = value;
        _high = value + (value < _bound ? _abs_lim : _rel_lim * value);
    }
    if (value > _max) {
        _max = value;
    }

    _sum += value;
    _sum_sq = std::fma(value, value, _sum_sq);
}

void GvcfBlock::reset()
{
    _min = INFINITY;
    _max = 0.0;
    _high = 0.0;
    _sum = 0.0;
    _sum_sq = 0.0;
}

} // namepace xatlas
