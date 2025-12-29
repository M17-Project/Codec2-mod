#include "codec2_internal.h"
#include "interp.h"
#include <math.h>

void interp_Wo(
    model_t *interp, /* interpolated model params */
    model_t *prev,   /* previous frames model params */
    model_t *next    /* next frames model params */
)
{
    const float weight = 0.5f;

    /* trap corner case where voicing est is probably wrong */
    if (interp->voiced && !prev->voiced && !next->voiced)
    {
        interp->voiced = 0;
    }

    /* Wo depends on voicing of this and adjacent frames */
    if (interp->voiced)
    {
        if (prev->voiced && next->voiced)
            interp->Wo = (1.0 - weight) * prev->Wo + weight * next->Wo;
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
    float *interp, /* interpolated LSPs */
    float *prev,   /* prev LSPs */
    float *next    /* next LSPs */
)
{
    const float weight = 0.5f;

    for (int i = 0; i < LPC_ORD; i++)
        interp[i] = (1.0f - weight) * prev[i] + weight * next[i];
}
