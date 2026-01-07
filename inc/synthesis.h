#ifndef CODEC2_MOD_SYNTHESIS_H
#define CODEC2_MOD_SYNTHESIS_H

#include <stdint.h>
#include "codec2_internal.h"

void synthesise_one_frame(
    codec2_t *c2,
    int16_t *speech,
    model_t *model,
    const complex_t *Aw,
    float gain
);

void synthesis_init(codec2_t *c2);

#endif /* CODEC2_MOD_SYNTHESIS_H */