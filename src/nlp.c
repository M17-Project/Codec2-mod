#include "codec2_internal.h"
#include "nlp.h"
#include <string.h>

/* original taps, but re-ordered in a polyphase filter bank */
const float nlp_fir_poly[NLP_NPHASE * NLP_TAPS_PP] = {
    // phase 0
    -0.00108181240f,
    0.00200298490f,
    0.00080284511f,
    -0.02092061000f,
    0.08630594400f,
    0.13674206000f,
    0.00322047750f,
    -0.01170581000f,
    0.00514494150f,
    -0.00092768838f,

    // phase 1
    -0.00110083440f,
    0.00370585090f,
    -0.00482046100f,
    -0.01280883100f,
    0.11480192000f,
    0.11480192000f,
    -0.01280883100f,
    -0.00482046100f,
    0.00370585090f,
    -0.00110083440f,

    // phase 2
    -0.00092768838f,
    0.00514494150f,
    -0.01170581000f,
    0.00322047750f,
    0.13674206000f,
    0.08630594400f,
    -0.02092061000f,
    0.00080284511f,
    0.00200298490f,
    -0.00108181240f,

    // phase 3
    -0.00042289438f,
    0.00559246660f,
    -0.01819927500f,
    0.02668381100f,
    0.14867556000f,
    0.05552062400f,
    -0.02206528200f,
    0.00430367540f,
    0.00055034190f,
    0.00000000000f,

    // phase 4
    0.00055034190f,
    0.00430367540f,
    -0.02206528200f,
    0.05552062400f,
    0.14867556000f,
    0.02668381100f,
    -0.01819927500f,
    0.00559246660f,
    -0.00042289438f,
    0.00000000000f,
};

static float post_process_sub_multiples(complex_t *Fw, float gmax, int gmax_bin, float *prev_f0)
{
    int min_bin, cmax_bin;
    int mult;
    float thresh, best_f0;
    int b, bmin, bmax, lmax_bin;
    float lmax;
    int prev_f0_bin;

    /* post process estimate by searching submultiples */
    mult = 2;
    min_bin = PE_FFT_SIZE * DEC / P_MAX;
    cmax_bin = gmax_bin;
    prev_f0_bin = *prev_f0 * (PE_FFT_SIZE * DEC) / SAMP_RATE;

    while (gmax_bin / mult >= min_bin)
    {
        b = gmax_bin / mult; /* determine search interval */
        bmin = 0.8 * b;
        bmax = 1.2 * b;
        if (bmin < min_bin)
            bmin = min_bin;

        /* lower threshold to favour previous frames pitch estimate,
            this is a form of pitch tracking */
        if ((prev_f0_bin > bmin) && (prev_f0_bin < bmax))
            thresh = CNLP * 0.5 * gmax;
        else
            thresh = CNLP * gmax;

        lmax = 0;
        lmax_bin = bmin;
        for (b = bmin; b <= bmax; b++) /* look for maximum in interval */
        {
            if (Fw[b].r > lmax)
            {
                lmax = Fw[b].r;
                lmax_bin = b;
            }
        }

        if (lmax > thresh)
        {
            if (lmax_bin > 0 && lmax_bin < (PE_FFT_SIZE / 2) &&
                lmax > Fw[lmax_bin - 1].r && lmax > Fw[lmax_bin + 1].r)
            {
                cmax_bin = lmax_bin;
            }
        }

        mult++;
    }

    best_f0 = (float)cmax_bin * SAMP_RATE / (PE_FFT_SIZE * DEC);

    return best_f0;
}

float nlp(
    nlp_t *restrict nlp,
    const float *restrict Sn, /* input speech vector */
    float *restrict pitch,    /* estimated pitch period in samples at current Fs    */
    float *restrict prev_f0   /* previous pitch f0 in Hz, memory for pitch tracking */
)
{
    float notch;                      /* current notch filter output          */
    complex_t *restrict Fw = nlp->Fw; /* DFT of squared signal (input/output) */
    float gmax;
    int gmax_bin;
    float best_f0;

    const int m = M_PITCH;
    const int n = N_SAMP;
    static const int d = N_SAMP / DEC;
    static const int START_POS = M_PITCH - N_SAMP;

    /* Square, notch filter at DC, and LP filter vector */
    /* Square latest input samples */
    for (int i = START_POS; i < m; i++)
    {
        nlp->sq[i] = Sn[i] * Sn[i];
    }

    for (int i = START_POS; i < m; i++)
    { /* notch filter at DC */
        notch = nlp->sq[i] - nlp->mem_x;
        notch += COEFF * nlp->mem_y;
        nlp->mem_x = nlp->sq[i];
        nlp->mem_y = notch;
        nlp->sq[i] = notch + 1.0; /* With 0 input vectors to codec,
                                     kiss_fft() would take a long
                                     time to execute when running in
                                     real time.  Problem was traced
                                     to kiss_fft function call in
                                     this function. Adding this small
                                     constant fixed problem.  Not
                                     exactly sure why. */
    }

    /* decimating polyphase FIR filter */
    /* push left */
    memmove(nlp->sq_fir, nlp->sq_fir + d, (NDEC - d) * sizeof(float));

    for (int i = 0; i < d; i++)
    {
        int idx = START_POS + NLP_NPHASE * i;
        float acc = 0.0f;

        for (int p = 0; p < NLP_NPHASE; p++)
        {
            const float *restrict h = &nlp_fir_poly[p * NLP_TAPS_PP];
            int base = idx - p;

            acc += h[0] * nlp->sq[base - 0 * NLP_NPHASE];
            acc += h[1] * nlp->sq[base - 1 * NLP_NPHASE];
            acc += h[2] * nlp->sq[base - 2 * NLP_NPHASE];
            acc += h[3] * nlp->sq[base - 3 * NLP_NPHASE];
            acc += h[4] * nlp->sq[base - 4 * NLP_NPHASE];
            acc += h[5] * nlp->sq[base - 5 * NLP_NPHASE];
            acc += h[6] * nlp->sq[base - 6 * NLP_NPHASE];
            acc += h[7] * nlp->sq[base - 7 * NLP_NPHASE];
            acc += h[8] * nlp->sq[base - 8 * NLP_NPHASE];
            acc += h[9] * nlp->sq[base - 9 * NLP_NPHASE];
        }

        nlp->sq_fir[NDEC - d + i] = acc; /* overwrite with d=16 fresh samples */
    }

    /* Decimate and DFT */
    memset(nlp->fftr_buff, 0, sizeof(nlp->fftr_buff));

    for (int i = 0; i < NDEC; i++)
    {
        nlp->fftr_buff[i] = nlp->sq_fir[i] * nlp->w[i];
    }

    kiss_fftr(nlp->fftr_cfg, nlp->fftr_buff, Fw);

    for (int i = 0; i < PE_FFT_SIZE / 2 + 1; i++)
    {
        Fw[i].r = Fw[i].r * Fw[i].r + Fw[i].i * Fw[i].i;
    }

    /* todo: express everything in f0, as pitch in samples is dep on Fs */
    int pmin = P_MIN;
    int pmax = P_MAX;

    /* find global peak */
    gmax = 0.0f;
    gmax_bin = PE_FFT_SIZE * DEC / pmax;
    int lo = PE_FFT_SIZE * DEC / pmax;
    int hi = PE_FFT_SIZE * DEC / pmin;
    if (hi > PE_FFT_SIZE / 2)
        hi = PE_FFT_SIZE / 2;

    for (int i = lo; i <= hi; i++)
    {
        if (Fw[i].r > gmax)
        {
            gmax = Fw[i].r;
            gmax_bin = i;
        }
    }

    best_f0 = post_process_sub_multiples(Fw, gmax, gmax_bin, prev_f0);

    /* Shift samples in buffer to make room for new samples */
    for (int i = 0; i < START_POS; i++)
        nlp->sq[i] = nlp->sq[i + n];

    /* return pitch period in samples and F0 estimate */
    *pitch = (float)SAMP_RATE / best_f0;

    *prev_f0 = best_f0;

    return (best_f0);
}

void nlp_init(nlp_t *nlp)
{
    memset(nlp->sq_fir, 0, sizeof(nlp->sq_fir));

    for (int i = 0; i < NDEC; i++)
    {
        nlp->w[i] = 0.5 - 0.5 * cosf(TWO_PI * i / (NDEC - 1));
    }

    memset(nlp->sq, 0, sizeof(nlp->sq));

    nlp->mem_x = 0.0;
    nlp->mem_y = 0.0;
}
