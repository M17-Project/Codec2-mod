#ifndef CODEC2_MOD_NLP_H
#define CODEC2_MOD_NLP_H

#include "codec2_internal.h"

float nlp(
    nlp_t *restrict nlp,
    const float *restrict Sn, /* input speech vector */
    float *restrict pitch,    /* estimated pitch period in samples at current Fs    */
    float *restrict prev_f0   /* previous pitch f0 in Hz, memory for pitch tracking */
);

void nlp_init(nlp_t *nlp);

#endif /* CODEC2_MOD_NLP_H */
