#include "codec2_internal.h"
#include "nlp.h"
#include <string.h>

const float nlp_fir[NLP_NTAP] = {
    -1.0818124e-03, -1.1008344e-03, -9.2768838e-04, -4.2289438e-04,
    5.5034190e-04, 2.0029849e-03, 3.7058509e-03, 5.1449415e-03,
    5.5924666e-03, 4.3036754e-03, 8.0284511e-04, -4.8204610e-03,
    -1.1705810e-02, -1.8199275e-02, -2.2065282e-02, -2.0920610e-02,
    -1.2808831e-02, 3.2204775e-03, 2.6683811e-02, 5.5520624e-02,
    8.6305944e-02, 1.1480192e-01, 1.3674206e-01, 1.4867556e-01,
    1.4867556e-01, 1.3674206e-01, 1.1480192e-01, 8.6305944e-02,
    5.5520624e-02, 2.6683811e-02, 3.2204775e-03, -1.2808831e-02,
    -2.0920610e-02, -2.2065282e-02, -1.8199275e-02, -1.1705810e-02,
    -4.8204610e-03, 8.0284511e-04, 4.3036754e-03, 5.5924666e-03,
    5.1449415e-03, 3.7058509e-03, 2.0029849e-03, 5.5034190e-04,
    -4.2289438e-04, -9.2768838e-04, -1.1008344e-03, -1.0818124e-03};

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
    int n,                    /* frames shift (no. new samples in Sn[])             */
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

    /* Square, notch filter at DC, and LP filter vector */
    /* Square latest input samples */
    for (int i = m - n; i < m; i++)
    {
        nlp->sq[i] = Sn[i] * Sn[i];
    }

    for (int i = m - n; i < m; i++)
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

    /* FIR filter vector */
    float *restrict mem = nlp->mem_fir;
    const float *restrict fir = nlp_fir;
    int pos = nlp->mem_pos;

    for (int i = m - n; i < m; i++)
    {
        /* write new sample twice */
        mem[pos] = nlp->sq[i];
        mem[pos + NLP_NTAP] = nlp->sq[i];

        /* FIR dot product: identical order to original */
        float acc = 0.0f;
        float *restrict x = &mem[pos + 1];

        /* NLP_NTAP=48 taps, unrolled by 4 */
        for (int j = 0; j < NLP_NTAP; j += 4)
        {
            acc += x[j + 0] * fir[j + 0];
            acc += x[j + 1] * fir[j + 1];
            acc += x[j + 2] * fir[j + 2];
            acc += x[j + 3] * fir[j + 3];
        }

        nlp->sq[i] = acc;

        /* advance pointer, wrap manually */
        pos++;
        if (pos == NLP_NTAP)
            pos = 0;
    }

    nlp->mem_pos = pos;

    /* Decimate and DFT */
    memset(nlp->fftr_buff, 0, sizeof(nlp->fftr_buff));

    for (int i = 0; i < NDEC; i++)
    {
        nlp->fftr_buff[i] = nlp->sq[i * DEC] * nlp->w[i];
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
    for (int i = 0; i < m - n; i++)
        nlp->sq[i] = nlp->sq[i + n];

    /* return pitch period in samples and F0 estimate */
    *pitch = (float)SAMP_RATE / best_f0;

    *prev_f0 = best_f0;

    return (best_f0);
}

void nlp_init(nlp_t *nlp)
{
	/* mem_fir must be cleared whenever mem_pos is reset */
	nlp->mem_pos = 0;
	memset(nlp->mem_fir, 0, sizeof(nlp->mem_fir));

	for (int i = 0; i < NDEC; i++)
	{
		nlp->w[i] = 0.5 - 0.5 * cosf(TWO_PI * i / (NDEC - 1));
	}

	memset(nlp->sq, 0, sizeof(nlp->sq));

	nlp->mem_x = 0.0;
	nlp->mem_y = 0.0;
}
