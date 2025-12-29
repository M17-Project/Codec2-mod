#ifndef INTERP_H
#define INTERP_H

void interp_Wo(model_t *interp, model_t *prev, model_t *next);
float interp_energy(float prev_e, float next_e);
void interpolate_lsp(float *interp, float *prev, float *next);

#endif /* INTERP_H */
