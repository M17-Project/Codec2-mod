#include "codec2_internal.h"
#include "util.h"
#include "lpc.h"

/* ~15Hz bandwidth expansion */
const float bw_gamma[LPC_ORD + 1] = {
    1.000000000000000,
    0.994000000000000,
    0.988036000000000,
    0.982107784000000,
    0.976215137296000,
    0.970357846472224,
    0.964535699393391,
    0.958748485197030,
    0.952995994285848,
    0.947278018320133,
    0.941594350210212};

static void autocorrelate(
    const float *restrict Sn, /* frame of Nsam windowed speech samples */
    float *restrict Rn        /* array of P+1 autocorrelation coefficients */
)
{
    for (int j = 0; j < LPC_ORD + 1; j++)
    {
        Rn[j] = 0.0;
        for (int i = 0; i < M_PITCH - j; i++)
            Rn[j] += Sn[i] * Sn[i + j];
    }
}

static void levinson_durbin(
    const float *restrict R, /* order+1 autocorrelation coeff */
    float *restrict lpcs     /* order+1 LPC's */
)
{
    float a[LPC_ORD + 1];
    float a_prev[LPC_ORD + 1];
    float e, k, sum;

    /* a[0] = 1 by definition */
    a[0] = 1.0f;
    e = R[0];

    for (int i = 1; i <= LPC_ORD; i++)
    {
        sum = 0.0f;
        for (int j = 1; j <= i - 1; j++)
            sum += a_prev[j] * R[i - j];

        k = -(R[i] + sum) / e;

        if (fabsf(k) > 1.0f)
            k = 0.0f;

        a[i] = k;

        for (int j = 1; j <= i - 1; j++)
            a[j] = a_prev[j] + k * a_prev[i - j];

        e *= (1.0f - k * k);

        /* copy for next iteration */
        for (int j = 1; j <= i; j++)
            a_prev[j] = a[j];
    }

    /* output LPCs */
    lpcs[0] = 1.0f;
    for (int i = 1; i <= LPC_ORD; i++)
        lpcs[i] = a[i];
}

/* NOTE: this function uses fully unrolled loop for LPC_ORD=10 */
static inline float cheb_poly_eva(
    const float *restrict coef, /* coefficients of the polynomial to be evaluated */
    float x                     /* the point where polynomial is to be evaluated */
)
{
    // N = 5 (LPC_ORD/2)
    float T0 = 1.0f;
    float T1 = x;

    float sum = coef[5] * T0 + coef[4] * T1;

    const float two_x = 2.0f * x;

    // i = 2
    float T2 = two_x * T1 - T0;
    sum += coef[3] * T2;

    // i = 3
    float T3 = two_x * T2 - T1;
    sum += coef[2] * T3;

    // i = 4
    float T4 = two_x * T3 - T2;
    sum += coef[1] * T4;

    // i = 5
    float T5 = two_x * T4 - T3;
    sum += coef[0] * T5;

    return sum;
}

/* NOTE: This function uses 6 bisections */
static int lpc_to_lsp(
    const float *a, /* LPC coefficients */
    float *freq     /* LSP frequencies in radians */
)
{
    float psuml, psumr, psumm, xl, xr, xm;
    float *pt;     /* ptr used for cheb_poly_eval(), whether P' or Q' */
    int roots = 0; /* number of roots found */
    float Q[LPC_ORD + 1];
    float P[LPC_ORD + 1];

    int m = LPC_ORD / 2; /* order of P'(z) & Q'(z) polynimials 	*/

    P[0] = Q[0] = 1.0f;
    for (int i = 1; i <= m; i++)
    {
        P[i] = a[i] + a[LPC_ORD + 1 - i] - P[i - 1];
        Q[i] = a[i] - a[LPC_ORD + 1 - i] + Q[i - 1];
    }

    for (int i = 0; i < m; i++)
    {
        P[i] *= 2.0f;
        Q[i] *= 2.0f;
    }

    /* Search for a zero in P'(z) polynomial first and then alternate to Q'(z).
    Keep alternating between the two polynomials as each zero is found 	*/
    xl = 1.0; /* start at point xl = 1 		*/
    const float delta = LSP_DELTA1;

    for (int j = 0; j < LPC_ORD; j++)
    {
        pt = (j & 1) ? Q : P; /* determines whether P' or Q' is eval. */

        xr = xl;
        psuml = cheb_poly_eva(pt, xl); /* evals poly. at xl 	*/

        while (xr >= -1.0f)
        {
            xr = xl - delta;               /* interval spacing 	*/
            psumr = cheb_poly_eva(pt, xr); /* poly(xl-delta_x) 	*/

            /* if no sign change increment xr and re-evaluate
               poly(xr). Repeat til sign change.  if a sign change has
               occurred the interval is bisected and then checked again
               for a sign change which determines in which interval the
               zero lies in.  If there is no sign change between poly(xm)
               and poly(xl) set interval between xm and xr else set
               interval between xl and xr and repeat till root is located
               within the specified limits  */
            if ((psumr <= 0.0f && psuml >= 0.0f) || (psumr >= 0.0f && psuml <= 0.0f)) // avoid one float multiplication
            {
                roots++;

                // manually unrolled for nb=5 (thus 6x)
                xm = 0.5f * (xl + xr);
                psumm = cheb_poly_eva(pt, xm);
                (psumm * psuml > 0.f) ? (xl = xm, psuml = psumm) : (xr = xm);

                xm = 0.5f * (xl + xr);
                psumm = cheb_poly_eva(pt, xm);
                (psumm * psuml > 0.f) ? (xl = xm, psuml = psumm) : (xr = xm);

                xm = 0.5f * (xl + xr);
                psumm = cheb_poly_eva(pt, xm);
                (psumm * psuml > 0.f) ? (xl = xm, psuml = psumm) : (xr = xm);

                xm = 0.5f * (xl + xr);
                psumm = cheb_poly_eva(pt, xm);
                (psumm * psuml > 0.f) ? (xl = xm, psuml = psumm) : (xr = xm);

                xm = 0.5f * (xl + xr);
                psumm = cheb_poly_eva(pt, xm);
                (psumm * psuml > 0.f) ? (xl = xm, psuml = psumm) : (xr = xm);

                xm = 0.5f * (xl + xr);
                psumm = cheb_poly_eva(pt, xm);
                (psumm * psuml > 0.f) ? (xl = xm, psuml = psumm) : (xr = xm);

                /* once zero is found, reset initial interval to xr 	*/
                freq[j] = (xm);
                xl = xm;
                break;
            }
            else
            {
                psuml = psumr;
                xl = xr;
            }
        }
    }

    /* convert from x domain to radians */
    for (int i = 0; i < LPC_ORD; i++)
    {
        freq[i] = fast_acosf(freq[i]);
    }

    return (roots);
}

void lsp_to_lpc(
    const float *restrict lsp, /* array of LSP frequencies in radians */
    float *restrict ak         /* array of LPC coefficients */
)
{
    float xout1, xout2, xin1, xin2;
    float *n1, *n2, *n3, *n4 = 0;
    float freq[LPC_ORD];
    float Wp[2 * LPC_ORD + 2];

    /* convert from radians to the x=cos(w) domain */
    for (int i = 0; i < LPC_ORD; i++)
        freq[i] = fast_cosf(lsp[i]);

    /* initialise contents of array */
    memset(Wp, 0, sizeof(Wp));

    xin1 = 1.0;
    xin2 = 1.0;

    /* reconstruct P(z) and Q(z) by cascading second order polynomials
      in form 1 - 2xz(-1) +z(-2), where x is the LSP coefficient */
    for (int j = 0; j <= LPC_ORD; j++)
    {
        for (int i = 0; i < (LPC_ORD / 2); i++)
        {
            n1 = Wp + (i * 4);
            n2 = n1 + 1;
            n3 = n2 + 1;
            n4 = n3 + 1;
            xout1 = xin1 - 2 * (freq[2 * i]) * *n1 + *n2;
            xout2 = xin2 - 2 * (freq[2 * i + 1]) * *n3 + *n4;
            *n2 = *n1;
            *n4 = *n3;
            *n1 = xin1;
            *n3 = xin2;
            xin1 = xout1;
            xin2 = xout2;
        }

        xout1 = xin1 + *(n4 + 1);
        xout2 = xin2 - *(n4 + 2);

        ak[j] = (xout1 + xout2) * 0.5;

        *(n4 + 1) = xin1;
        *(n4 + 2) = xin2;

        xin1 = 0.0;
        xin2 = 0.0;
    }
}

void aks_to_mag2(codec2_t *c2,
                        const float *ak, /* LPCs */
                        model_t *model,  /* sinusoidal model parameters for this frame */
                        float E,         /* energy term */
                        complex_t *Aw,   /* output power spectrum */
                        float *A2)
{
    int am, bm; /* limits of current band */

    /* FFT of A(z) */
    float *a = (float *)c2->fft_buffer;
    memset(a, 0, FFT_ENC * sizeof(float));

    for (int i = 0; i <= LPC_ORD; i++)
        a[i] = ak[i];

    kiss_fftr(c2->fftr_fwd_cfg, a, Aw);

    for (int i = 0; i < FFT_ENC / 2; i++)
    {
        A2[i] = Aw[i].r * Aw[i].r + Aw[i].i * Aw[i].i + 1e-6f;
    }

    /* build ak_gamma */
    float *ag = (float *)c2->fft_buffer;
    memset(ag, 0, FFT_ENC * sizeof(float));

    float g = LPCPF_GAMMA;
    ag[0] = ak[0];
    for (int i = 1; i <= LPC_ORD; i++)
    {
        ag[i] = ak[i] * g;
        g *= LPCPF_GAMMA;
    }

    /* FFT of A_gamma */
    complex_t *Awg = c2->fft_buffer; /* reuse FFT scratch */
    kiss_fftr(c2->fftr_fwd_cfg, ag, Awg);

    /* reuse Awg storage for A2g */
    float *A2g = (float *)Awg;

    for (int i = 0; i < FFT_ENC / 2; i++)
    {
        A2g[i] = Awg[i].r * Awg[i].r + Awg[i].i * Awg[i].i + 1e-6f;
    }

    /* compute normalization gain */
    float e_before = 1e-12f;
    float e_after = 1e-12f;

    for (int i = 0; i < FFT_ENC / 2; i++)
    {
        float invA2 = 1.0f / A2[i];
        float R = sqrtf(A2g[i] * invA2);

        e_before += invA2;
        e_after += invA2 * expf(LPCPF_TWO_BETA * logf(R + 1e-5f));
    }

    float gain = E * e_before / e_after;

    /* Determine magnitudes */
    for (int m = 1; m <= model->L; m++)
    {
        am = (int)((m - 0.5f) * model->Wo / FFT_R + 0.5f);
        bm = (int)((m + 0.5f) * model->Wo / FFT_R + 0.5f);

        if (bm > FFT_ENC / 2)
            bm = FFT_ENC / 2;

        float Em = 0.0f;

        for (int i = am; i < bm; i++)
        {
            /* R(w) = |A_gamma| / |A| */
            float R = sqrtf(A2g[i] / A2[i]);

            /* Pw contribution */
            float Pw_i = expf(LPCPF_TWO_BETA * logf(R + 1e-5f)) / A2[i];

            /* boost low frequencies a bit */
            float freq = i * (SAMP_RATE * 0.5f / (FFT_ENC / 2));
            if (freq < 1000.0f)
                Pw_i *= 1.96f;

            Em += Pw_i;
        }

        /* apply LPC energy */
        Em *= gain;

        model->A[m] = sqrtf(Em);
    }
}

/* Apply first harmonic LPC correction at decoder.
   This helps improve low pitch males after LPC modelling. */
void apply_lpc_correction(model_t *model)
{
    if (model->Wo < (M_PI * 150.0f / 4000.0f))
    {
        model->A[1] *= 0.032f;
    }
}

float speech_to_uq_lsps(
    codec2_t *c2,
    float *restrict lsp,
    float *restrict ak,
    float *restrict energy,
    const float *restrict Sn,
    const float *restrict w /* window */
)
{
    int roots;
    float R[LPC_ORD + 1];
    float e, E;

    float *Wn = (float *)c2->fft_buffer;

    e = 0.0;
    for (int i = 0; i < M_PITCH; i++)
    {
        Wn[i] = Sn[i] * w[i];
        e += Wn[i] * Wn[i];
    }

    /* trap 0 energy case as LPC analysis will fail */
    if (e < LPC_ENERGY_FLOOR)
    {
        for (int i = 0; i < LPC_ORD; i++)
            lsp[i] = (M_PI / LPC_ORD) * (float)i;

        memset(ak, 0, (LPC_ORD + 1) * sizeof(float));
        ak[0] = 1.0f;
        *energy = 0.0f;
        return 0.0f;
    }

    autocorrelate(Wn, R);
    levinson_durbin(R, ak);

    E = 0.0f;
    for (int i = 0; i <= LPC_ORD; i++)
        E += ak[i] * R[i];

    /* 15 Hz BW expansion as I can't hear the difference and it may help
       help occasional fails in the LSP root finding.  Important to do this
       after energy calculation to avoid -ve energy values.
    */
    for (int i = 0; i <= LPC_ORD; i++)
    {
        ak[i] *= bw_gamma[i];
    }

    roots = lpc_to_lsp(ak, lsp); // hardcoded to 5 bisections
    if (roots != LPC_ORD)
    {
        /* if root finding fails use some benign LSP values instead */
        for (int i = 0; i < LPC_ORD; i++)
            lsp[i] = (M_PI / LPC_ORD) * (float)i;
    }

    *energy = E;

    return (E >= 0.0f) ? 0 : 1;
}
