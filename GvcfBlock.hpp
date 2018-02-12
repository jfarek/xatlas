#ifndef _XATLAS_BLOCKSTATS_H
#define _XATLAS_BLOCKSTATS_H

#include <cstdint>

class GvcfBlock
{
  private:
    double _high;    // upper bound
    double _abs_lim; // absolute limit
    double _rel_lim; // relative limit coefficient
    double _sum;     // sum of values
    double _sum_sq;  // sum of squared values
    double _bound;   // rel/abs limit boundary

    GvcfBlock();

  public:
    double _min; // minimum value
    double _max; // maximum value

    GvcfBlock(double block_abs_lim, double block_rel_lim);
    GvcfBlock(const GvcfBlock &);
    ~GvcfBlock();
    bool break_block(double value);
    void add_value(double value);
    void reset();
    double get_std_dev(uint32_t k);
    double get_mean(uint32_t k);
};

#endif /* _XATLAS_BLOCKSTATS_H */
