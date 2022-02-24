#include <cmath>

#include "xatlas/Logit.h"

namespace xatlas {

// \frac{1}{1/e^{-x}}

double snp_logit(snp_logit_params_s *params, snp_logit_values_s *values)
{
    return 1.0 / (1.0 + std::exp(-std::fma(params->ratio_score, values->ratio_score,
                                 std::fma(params->base_qual,    values->base_qual,
                                 std::fma(params->mean_avnqs,   values->mean_avnqs,
                                 std::fma(params->rel_pos,      values->rel_pos,
                                 std::fma(params->titv,         values->titv,
                                 params->intercept)))))));
}

double indel_logit(indel_logit_params_s *params, indel_logit_values_s *values)
{
    return 1.0 / (1.0 + std::exp(-std::fma(params->ratio_score,  values->ratio_score,
                                 std::fma(params->strand_dir,    values->strand_dir,
                                 std::fma(params->mean_avnqs,    values->mean_avnqs,
                                 std::fma(params->seq_entropy,   values->seq_entropy,
                                 std::fma(params->mean_var_rate, values->mean_var_rate,
                                 params->intercept)))))));
}

} // namepace xatlas
