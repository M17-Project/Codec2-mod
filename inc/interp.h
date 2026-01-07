#ifndef CODEC2_MOD_INTERP_H
#define CODEC2_MOD_INTERP_H

#include "codec2_internal.h"

void interp_Wo(model_t *interp, const model_t *prev, const model_t *next);
float interp_energy(float prev_e, float next_e);
void interpolate_lsp(float *interp, const float *prev, const float *next);

#endif /* CODEC2_MOD_INTERP_H */
