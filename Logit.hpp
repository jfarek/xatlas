#ifndef _XATLAS_LOGIT_H
#define _XATLAS_LOGIT_H

/* logit model parameters */

struct snp_logit_params {
    double intercept;
    double ratio_score;
    double base_qual;
    double mean_avnqs;
    double rel_pos;
    double titv;
};
typedef struct snp_logit_params snp_logit_params_s;

struct indel_logit_params {
    double intercept;
    double ratio_score;
    double strand_dir;
    double mean_avnqs;
    double seq_entropy;
    double mean_var_rate;
};
typedef struct indel_logit_params indel_logit_params_s;

struct logit_params {
    struct snp_logit_params *snp;
    struct indel_logit_params *indel;
};
typedef struct logit_params logit_params_s;

/* values */

struct snp_logit_values {
    double ratio_score;
    double base_qual;
    double mean_avnqs;
    double rel_pos;
    double titv;
};
typedef struct snp_logit_values snp_logit_values_s;

struct indel_logit_values {
    double ratio_score;
    double strand_dir;
    double mean_avnqs;
    double seq_entropy;
    double mean_var_rate;
};
typedef struct indel_logit_values indel_logit_values_s;

double snp_logit(snp_logit_params_s *params, snp_logit_values_s *values);
double indel_logit(indel_logit_params_s *params, indel_logit_values_s *values);

#endif /* _XATLAS_LOGIT_H */
