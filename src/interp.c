#include "codec2_internal.h"
#include "interp.h"
#include <math.h>

void interp_Wo(
    model_t *interp,              /* interpolated model params */
    const model_t *restrict prev, /* previous frames model params */
    const model_t *restrict next  /* next frames model params */
)
{
    /* trap corner case where voicing est is probably wrong */
    if (interp->voiced && !prev->voiced && !next->voiced)
    {
        interp->voiced = 0;
    }

    /* Wo depends on voicing of this and adjacent frames */
    if (interp->voiced)
    {
        if (prev->voiced && next->voiced)
            interp->Wo = prev->Wo + 0.5f * (next->Wo - prev->Wo); // simple average, but disguised
        if (!prev->voiced && next->voiced)
            interp->Wo = next->Wo;
        if (prev->voiced && !next->voiced)
            interp->Wo = prev->Wo;
    }
    else
    {
        interp->Wo = W0_MIN;
    }

    interp->L = M_PI / interp->Wo;
}

float interp_energy(
    float prev_e, /* previous energy */
    float next_e  /* next energy */
)
{
    return sqrtf(prev_e * next_e);
}

void interpolate_lsp(
    float *interp,              /* interpolated LSPs */
    const float *restrict prev, /* prev LSPs */
    const float *restrict next  /* next LSPs */
)
{
    for (int i = 0; i < LPC_ORD; i++)
        interp[i] = prev[i] + 0.5f * (next[i] - prev[i]); // simple average, disguised
}
