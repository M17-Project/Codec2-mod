#ifndef CODEC2_MOD_LPC_H
#define CODEC2_MOD_LPC_H

void aks_to_mag2(
    codec2_t *c2,
    const float *ak,
    model_t *model,
    float E,
    complex_t *Aw,
    float *A2);

float speech_to_uq_lsps(
    codec2_t *c2,
    float *restrict lsp,
    float *restrict ak,
    float *restrict energy,
    const float *restrict Sn,
    const float *restrict w);

void lsp_to_lpc(
    const float *restrict lsp,
    float *restrict ak);

void apply_lpc_correction(model_t *model);

#endif /* CODEC2_MOD_LPC_H */
