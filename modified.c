/*
  This work uses code written by David Rowe VK5DGR et al.
  https://github.com/drowe67/codec2
*/

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

#include "kiss_fft.h"
#include "kiss_fftr.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TWO_PI (2.0f * M_PI)

#ifndef POW10F
#define POW10F(x) powf(10.0f, (x))
#endif

#define N_S 0.01   /* internal proc frame length in secs   */
#define MAX_AMP 80 /* maximum number of harmonics          */

#define FFT_ENC 512	 /* size of FFT used for encoder         */
#define FFT_DEC 512	 /* size of FFT used in decoder          */
#define V_THRESH 6.0 /* voicing threshold in dB              */
#define LPC_ORD 10	 /* LPC order                            */

/* Pitch estimation defines */
#define P_MAX_S 0.0200f /* maximum pitch period in s            */

#define SAMP_RATE 8000 /* sample rate of this instance             */
#define N_SAMP 80	   /* number of samples per 10ms frame at Fs   */
#define M_PITCH 320	   /* pitch estimation window size in samples  */
#define P_MIN 20	   /* minimum pitch period in samples          */
#define P_MAX 160	   /* maximum pitch period in samples          */
#define W0_MIN ((2.0 * M_PI) / 160.0)
#define W0_MAX ((2.0 * M_PI) / 20.0)
#define NW 279 /* analysis window size in samples          */
#define TW 40  /* trapezoidal synthesis window overlap     */

/* NLP */
#define PMAX_M 320		/* maximum NLP analysis window size     */
#define COEFF 0.95		/* notch filter parameter               */
#define PE_FFT_SIZE 512 /* DFT size for pitch estimation        */
#define DEC 5			/* decimation factor                    */
#define NDEC (M_PITCH / DEC)
#define T 0.1		/* threshold for local minima candidate */
#define CNLP 0.3	/* post processor constant              */
#define NLP_NTAP 48 /* Decimation LPF order */

/* quantizers & LPC */
#define WO_BITS 7
#define WO_LEVELS (1 << WO_BITS)
#define E_BITS 5
#define E_LEVELS (1 << E_BITS)
#define E_MIN_DB -10.0
#define E_MAX_DB 40.0
#define LSPD_SCALAR_INDEXES 10

#define LPCPF_GAMMA 0.5
#define LPCPF_BETA 0.2
#define LPCPF_TWO_BETA (2.0 * LPCPF_BETA)
#define LSP_DELTA1 0.01 /* grid spacing for LSP root searches */
#define BG_THRESH 40.0	/* only consider low levels signals for bg_est */
#define BG_BETA 0.1		/* averaging filter constant                   */
#define BG_MARGIN 6.0

/* FFT stuff */
#define FFT_R (TWO_PI / FFT_ENC) /* radians/bin */
#define FFT_1_R (FFT_ENC / TWO_PI)

/* random number generator stuff */
#define CODEC2_RAND_MAX 32767

/* consts */
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

extern const float delta_lsp_cb[10][32];

typedef kiss_fft_cpx complex_t;

typedef struct model_t
{
	float Wo;				/* fundamental frequency estimate in radians  */
	int L;					/* number of harmonics                        */
	float A[MAX_AMP + 1];	/* amplitiude of each harmonic                */
	float phi[MAX_AMP + 1]; /* phase of each harmonic                     */
	int voiced;				/* non-zero if this frame is voiced           */
} model_t;

typedef struct nlp_t
{
	float w[PMAX_M / DEC];		 /* DFT window                   */
	float sq[PMAX_M];			 /* squared speech samples       */
	float mem_x, mem_y;			 /* memory for notch filter      */
	float mem_fir[NLP_NTAP * 2]; /* decimation FIR filter memory */
	int mem_pos;
	kiss_fftr_cfg fftr_cfg; /* kiss FFT config              */
	float fftr_buff[PE_FFT_SIZE];
	complex_t Fw[PE_FFT_SIZE / 2 + 1];
} nlp_t;

typedef struct codec2_t
{
	unsigned long next_rn; /* for the pseudorandom number geneartor     */

	float w[M_PITCH];	  /* [m_pitch] time domain hamming window      */
	float W[FFT_ENC];	  /* DFT of w[]                                */
	float Pn[2 * N_SAMP]; /* [2*n_samp] trapezoidal synthesis window   */
	float Sn[M_PITCH];	  /* [m_pitch] input speech                    */
	float hpf_states[2];  /* high pass filter states                   */
	nlp_t nlp;			  /* pitch predictor states                    */

	float Sn_[2 * N_SAMP];		  /* [2*n_samp] synthesised output speech      */
	float ex_phase;				  /* excitation model phase track              */
	float bg_est;				  /* background noise estimate for post filter */
	float prev_f0_enc;			  /* previous frame's f0    estimate           */
	model_t prev_model_dec;		  /* previous frame's model parameters         */
	float prev_lsps_dec[LPC_ORD]; /* previous frame's LSPs                     */
	float prev_e_dec;			  /* previous frame's LPC energy               */

	int lpc_pf;		/* LPC post filter on                        */
	int bass_boost; /* LPC post filter bass boost                */

	kiss_fft_cfg fft_fwd_cfg;	/* forward FFT config                        */
	kiss_fftr_cfg fftr_fwd_cfg; /* forward real FFT config                   */
	kiss_fftr_cfg fftr_inv_cfg; /* inverse FFT config                        */
	kiss_fft_cfg phase_fft_fwd_cfg;
	kiss_fft_cfg phase_fft_inv_cfg;
	kiss_fft_cpx fft_buffer[FFT_ENC];
} codec2_t;

void make_analysis_window(kiss_fft_cfg fft_fwd_cfg, float *w, float *W)
{
	complex_t wshift[FFT_ENC] = {0};

	float m = 0.0;
	for (int i = 0; i < M_PITCH / 2 - NW / 2; i++)
		w[i] = 0.0f;
	for (int i = M_PITCH / 2 - NW / 2, j = 0; i < M_PITCH / 2 + NW / 2; i++, j++)
	{
		w[i] = 0.5f - 0.5f * cosf(TWO_PI * j / (NW - 1));
		m += w[i] * w[i];
	}
	for (int i = M_PITCH / 2 + NW / 2; i < M_PITCH; i++)
		w[i] = 0.0f;

	m = 1.0f / sqrtf(m * FFT_ENC);
	for (int i = 0; i < M_PITCH; i++)
	{
		w[i] *= m;
	}

	for (int i = 0; i < NW / 2; i++)
		wshift[i].r = w[i + M_PITCH / 2];

	for (int i = FFT_ENC - NW / 2, j = M_PITCH / 2 - NW / 2; i < FFT_ENC; i++, j++)
		wshift[i].r = w[j];

	kiss_fft(fft_fwd_cfg, wshift, wshift);

	for (int i = 0; i < FFT_ENC / 2; i++)
	{
		W[i] = wshift[i + FFT_ENC / 2].r;
		W[i + FFT_ENC / 2] = wshift[i].r;
	}
}

void make_synthesis_window(float *Pn)
{
	int i;
	float win;

	/* Generate Parzen window in time domain */
	win = 0.0;
	for (i = 0; i < N_SAMP / 2 - TW; i++)
		Pn[i] = 0.0;
	win = 0.0;
	for (i = N_SAMP / 2 - TW; i < N_SAMP / 2 + TW; win += 1.0 / (2 * TW), i++)
		Pn[i] = win;
	for (i = N_SAMP / 2 + TW; i < 3 * N_SAMP / 2 - TW; i++)
		Pn[i] = 1.0;
	win = 1.0;
	for (i = 3 * N_SAMP / 2 - TW; i < 3 * N_SAMP / 2 + TW;
		 win -= 1.0f / (2.0f * TW), i++)
		Pn[i] = win;
	for (i = 3.0f * N_SAMP / 2.0f + TW; i < 2 * N_SAMP; i++)
		Pn[i] = 0.0;
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

	nlp->fftr_cfg = kiss_fftr_alloc(PE_FFT_SIZE, 0, NULL, NULL);
}

void codec2_init(codec2_t *c2)
{
	_Static_assert(PE_FFT_SIZE % 2 == 0, "PE_FFT_SIZE must be even");

	c2->next_rn = 1; // random number geterator - seed

	for (int i = 0; i < M_PITCH; i++)
		c2->Sn[i] = 1.0f;

	c2->hpf_states[0] = c2->hpf_states[1] = 0.0;

	memset(c2->Sn_, 0, sizeof(c2->Sn_));

	c2->fft_fwd_cfg = kiss_fft_alloc(FFT_ENC, 0, NULL, NULL);
	c2->fftr_fwd_cfg = kiss_fftr_alloc(FFT_ENC, 0, NULL, NULL);

	make_analysis_window(c2->fft_fwd_cfg, c2->w, c2->W);
	make_synthesis_window(c2->Pn);

	c2->fftr_inv_cfg = kiss_fftr_alloc(FFT_DEC, 1, NULL, NULL);
	c2->prev_f0_enc = 1.0f / P_MAX_S;
	c2->bg_est = 0.0;
	c2->ex_phase = 0.0;

	for (int l = 1; l <= MAX_AMP; l++)
		c2->prev_model_dec.A[l] = 0.0;

	c2->prev_model_dec.Wo = TWO_PI / P_MAX;
	c2->prev_model_dec.L = M_PI / c2->prev_model_dec.Wo;
	c2->prev_model_dec.voiced = 0;
	memset(c2->prev_model_dec.phi, 0, sizeof(c2->prev_model_dec.phi));
	c2->ex_phase = 0.0f;

	for (int i = 0; i < LPC_ORD; i++)
	{
		c2->prev_lsps_dec[i] = i * M_PI / (LPC_ORD + 1);
	}
	c2->prev_e_dec = 1;

	nlp_init(&c2->nlp);

	c2->lpc_pf = 1;
	c2->bass_boost = 1;
}

void dft_speech(kiss_fft_cfg fft_fwd_cfg, complex_t *Sw, float *Sn, float *w)
{
	memset(Sw, 0, FFT_ENC * sizeof(*Sw));

	/* Centre analysis window on time axis, we need to arrange input
	   to FFT this way to make FFT phases correct */
	/* move 2nd half to start of FFT input vector */
	for (int i = 0; i < NW / 2; i++)
		Sw[i].r = Sn[i + M_PITCH / 2] * w[i + M_PITCH / 2];

	/* move 1st half to end of FFT input vector */
	for (int i = 0; i < NW / 2; i++)
		Sw[FFT_ENC - NW / 2 + i].r =
			Sn[i + M_PITCH / 2 - NW / 2] * w[i + M_PITCH / 2 - NW / 2];

	kiss_fft(fft_fwd_cfg, Sw, Sw);
}

void estimate_amplitudes(model_t *model, complex_t *Sw, int est_phase)
{
	int am, bm; /* bounds of current harmonic */
	float den;	/* denominator of amplitude expression */

	for (int m = 1; m <= model->L; m++)
	{
		/* Estimate ampltude of harmonic */
		den = 0.0f;
		am = (int)((m - 0.5f) * model->Wo * FFT_1_R + 0.5f);
		bm = (int)((m + 0.5f) * model->Wo * FFT_1_R + 0.5f);

		for (int i = am; i < bm; i++)
		{
			den += Sw[i].r * Sw[i].r + Sw[i].i * Sw[i].i;
		}

		model->A[m] = sqrtf(den);

		if (est_phase)
		{
			int b = (int)(m * model->Wo / FFT_R + 0.5); /* DFT bin of centre of current harmonic */

			/* Estimate phase of harmonic, this is expensive in CPU for
			   embedded devicesso we make it an option */
			model->phi[m] = atan2f(Sw[b].i, Sw[b].r);
		}
	}
}

float post_process_sub_multiples(complex_t *Fw, float gmax, int gmax_bin, float *prev_f0)
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
			if ((lmax > Fw[lmax_bin - 1].r) && (lmax > Fw[lmax_bin + 1].r))
			{
				cmax_bin = lmax_bin;
			}
		}

		mult++;
	}

	best_f0 = (float)cmax_bin * SAMP_RATE / (PE_FFT_SIZE * DEC);

	return best_f0;
}

void hs_pitch_refinement(model_t *model, complex_t *Sw, float pmin, float pmax, float pstep)
{
	int m;	   /* loop variable */
	int b;	   /* bin for current harmonic centre */
	float E;   /* energy for current pitch*/
	float Wo;  /* current "test" fundamental freq. */
	float Wom; /* Wo that maximises E */
	float Em;  /* mamimum energy */
	float p;   /* current pitch */

	/* Initialisation */
	model->L = M_PI / model->Wo; /* use initial pitch est. for L */
	Wom = model->Wo;
	Em = 0.0;

	/* Determine harmonic sum for a range of Wo values */
	for (p = pmin; p <= pmax; p += pstep)
	{
		E = 0.0;
		Wo = TWO_PI / p;

		float bFloat = Wo * FFT_1_R;
		float currentBFloat = bFloat;

		/* Sum harmonic magnitudes */
		for (m = 1; m <= model->L; m++)
		{
			b = (int)(currentBFloat + 0.5);
			E += Sw[b].r * Sw[b].r + Sw[b].i * Sw[b].i;
			currentBFloat += bFloat;
		}
		/* Compare to see if this is a maximum */
		if (E > Em)
		{
			Em = E;
			Wom = Wo;
		}
	}

	model->Wo = Wom;
}

void two_stage_pitch_refinement(model_t *model, complex_t *Sw)
{
	float pmin, pmax, pstep; /* pitch refinement minimum, maximum and step */

	/* Coarse refinement */
	pmax = TWO_PI / model->Wo + 5;
	pmin = TWO_PI / model->Wo - 5;
	pstep = 1.0;
	hs_pitch_refinement(model, Sw, pmin, pmax, pstep);

	/* Fine refinement */
	pmax = TWO_PI / model->Wo + 1;
	pmin = TWO_PI / model->Wo - 1;
	pstep = 0.25;
	hs_pitch_refinement(model, Sw, pmin, pmax, pstep);

	/* Limit range */
	if (model->Wo < TWO_PI / P_MAX)
		model->Wo = TWO_PI / P_MAX;
	if (model->Wo > TWO_PI / P_MIN)
		model->Wo = TWO_PI / P_MIN;

	model->L = floorf(M_PI / model->Wo);

	/* trap occasional round off issues with floorf() */
	if (model->Wo * model->L >= 0.95 * M_PI)
	{
		model->L--;
	}
}

float est_voicing_mbe(model_t *model, complex_t *Sw, float *W)
{
	int l, al, bl, m; /* loop variables */
	complex_t Am;	  /* amplitude sample for this band */
	int offset;		  /* centers Hw[] about current harmonic */
	float den;		  /* denominator of Am expression */
	float error;	  /* accumulated error between original and synthesised */
	float Wo;
	float sig, snr;
	float elow, ehigh, eratio;
	float sixty;
	complex_t Ew;
	Ew.r = 0;
	Ew.i = 0;

	int l_1000hz = model->L * 1000.0 / (SAMP_RATE / 2);
	sig = 1E-4;
	for (l = 1; l <= l_1000hz; l++)
	{
		sig += model->A[l] * model->A[l];
	}

	Wo = model->Wo;
	error = 1E-4;

	/* Just test across the harmonics in the first 1000 Hz */
	for (l = 1; l <= l_1000hz; l++)
	{
		Am.r = 0.0;
		Am.i = 0.0;
		den = 0.0;
		al = ceilf((l - 0.5) * Wo * FFT_ENC / TWO_PI);
		bl = ceilf((l + 0.5) * Wo * FFT_ENC / TWO_PI);

		/* Estimate amplitude of harmonic assuming harmonic is totally voiced */
		offset = FFT_ENC / 2 - l * Wo * FFT_ENC / TWO_PI + 0.5;
		for (m = al; m < bl; m++)
		{
			Am.r += Sw[m].r * W[offset + m];
			Am.i += Sw[m].i * W[offset + m];
			den += W[offset + m] * W[offset + m];
		}

		Am.r = Am.r / den;
		Am.i = Am.i / den;

		/* Determine error between estimated harmonic and original */
		for (m = al; m < bl; m++)
		{
			Ew.r = Sw[m].r - Am.r * W[offset + m];
			Ew.i = Sw[m].i - Am.i * W[offset + m];
			error += Ew.r * Ew.r;
			error += Ew.i * Ew.i;
		}
	}

	snr = 10.0 * log10f(sig / error);
	if (snr > V_THRESH)
		model->voiced = 1;
	else
		model->voiced = 0;

	/* post processing, helps clean up some voicing errors ------------------*/
	/*
	   Determine the ratio of low frequency to high frequency energy,
	   voiced speech tends to be dominated by low frequency energy,
	   unvoiced by high frequency. This measure can be used to
	   determine if we have made any gross errors.
	*/
	int l_2000hz = model->L * 2000.0f / (SAMP_RATE / 2);
	int l_4000hz = model->L * 4000.0f / (SAMP_RATE / 2);
	elow = ehigh = 1E-4;
	for (l = 1; l <= l_2000hz; l++)
	{
		elow += model->A[l] * model->A[l];
	}
	for (l = l_2000hz; l <= l_4000hz; l++)
	{
		ehigh += model->A[l] * model->A[l];
	}
	eratio = 10.0 * log10f(elow / ehigh);

	/* Look for Type 1 errors, strongly V speech that has been
	   accidentally declared UV */
	if (model->voiced == 0)
		if (eratio > 10.0)
			model->voiced = 1;

	/* Look for Type 2 errors, strongly UV speech that has been
	   accidentally declared V */
	if (model->voiced == 1)
	{
		if (eratio < -10.0)
			model->voiced = 0;

		/* A common source of Type 2 errors is the pitch estimator
		   gives a low (50Hz) estimate for UV speech, which gives a
		   good match with noise due to the close harmoonic spacing.
		   These errors are much more common than people with 50Hz3
		   pitch, so we have just a small eratio threshold. */
		sixty = 60.0f * TWO_PI / SAMP_RATE;
		if ((eratio < -4.0f) && (model->Wo <= sixty))
			model->voiced = 0;
	}

	return snr;
}

float nlp(
	nlp_t *restrict nlp,
	float *restrict Sn,		/* input speech vector */
	int n,					/* frames shift (no. new samples in Sn[])             */
	float *restrict pitch,	/* estimated pitch period in samples at current Fs    */
	float *restrict prev_f0 /* previous pitch f0 in Hz, memory for pitch tracking */
)
{
	float notch;					  /* current notch filter output          */
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

void sample_phase(model_t *model, complex_t *H,
				  complex_t *A /* LPC analysis filter in freq domain */
)
{
	/* Sample phase at harmonics */
	for (int m = 1; m <= model->L; m++)
	{
		int b = (int)(m * model->Wo / FFT_R + 0.5);
		/* synth filter 1/A is opposite phase to analysis filter */
		// H[m] = cconj(A[b]);
		H[m].r = A[b].r;
		H[m].i = -A[b].i;
	}
}

int codec2_rand(unsigned long *prng_state)
{
	*prng_state = *prng_state * 1103515245 + 12345;
	return ((unsigned)(*prng_state / 65536) % 32768);
}

void postfilter(codec2_t *c2, model_t *model, float *bg_est)
{
	int m, uv;
	float e, thresh;

	/* determine average energy across spectrum */
	e = 1E-12;
	for (m = 1; m <= model->L; m++)
		e += model->A[m] * model->A[m];

	e = 10.0 * log10f(e / model->L);

	/* If beneath threshold, update bg estimate.  The idea
	   of the threshold is to prevent updating during high level
	   speech. */
	if ((e < BG_THRESH) && !model->voiced)
		*bg_est = *bg_est * (1.0 - BG_BETA) + e * BG_BETA;

	/* now mess with phases during voiced frames to make any harmonics
	   less then our background estimate unvoiced.
	*/
	uv = 0;
	thresh = POW10F((*bg_est + BG_MARGIN) / 20.0);
	if (model->voiced)
		for (m = 1; m <= model->L; m++)
			if (model->A[m] < thresh)
			{
				model->phi[m] = (TWO_PI / CODEC2_RAND_MAX) * (float)codec2_rand(&c2->next_rn);
				uv++;
			}
}

void synthesise(kiss_fftr_cfg fftr_inv_cfg,
				float *Sn_,		/* time domain synthesised signal              */
				model_t *model, /* ptr to model parameters for this frame      */
				float *Pn,		/* time domain Parzen window                   */
				int shift		/* flag used to handle transition frames       */
)
{
	int i, l, j, b;					/* loop variables */
	complex_t Sw_[FFT_DEC / 2 + 1]; /* DFT of synthesised signal */
	float sw_[FFT_DEC];				/* synthesised signal */

	if (shift)
	{
		/* Update memories */
		for (i = 0; i < N_SAMP - 1; i++)
		{
			Sn_[i] = Sn_[i + N_SAMP];
		}
		Sn_[N_SAMP - 1] = 0.0;
	}

	memset(Sw_, 0, sizeof(Sw_));

	/* Now set up frequency domain synthesised speech */
	for (l = 1; l <= model->L; l++)
	{
		b = (int)(l * model->Wo * FFT_DEC / TWO_PI + 0.5);
		if (b > ((FFT_DEC / 2) - 1))
		{
			b = (FFT_DEC / 2) - 1;
		}
		Sw_[b].r = model->A[l] * cosf(model->phi[l]);
		Sw_[b].i = model->A[l] * sinf(model->phi[l]);
	}

	/* Perform inverse DFT */
	kiss_fftri(fftr_inv_cfg, Sw_, sw_);

	/* Overlap add to previous samples */
	for (i = 0; i < N_SAMP - 1; i++)
	{
		Sn_[i] += sw_[FFT_DEC - N_SAMP + 1 + i] * Pn[i];
	}

	if (shift)
		for (i = N_SAMP - 1, j = 0; i < 2 * N_SAMP; i++, j++)
			Sn_[i] = sw_[j] * Pn[i];
	else
		for (i = N_SAMP - 1, j = 0; i < 2 * N_SAMP; i++, j++)
			Sn_[i] += sw_[j] * Pn[i];
}

void phase_synth_zero_order(
	codec2_t *c2,
	model_t *model,
	float *ex_phase, /* excitation phase of fundamental        */
	complex_t *H	 /* L synthesis filter freq domain samples */

)
{
	int m;
	float new_phi;
	complex_t Ex[MAX_AMP + 1]; /* excitation samples */
	complex_t A_[MAX_AMP + 1]; /* synthesised harmonic samples */

	/*
	   Update excitation fundamental phase track, this sets the position
	   of each pitch pulse during voiced speech.  After much experiment
	   I found that using just this frame's Wo improved quality for UV
	   sounds compared to interpolating two frames Wo like this:

	   ex_phase[0] += (*prev_Wo+model->Wo)*N_SAMP/2;
	*/
	ex_phase[0] += (model->Wo) * N_SAMP;
	ex_phase[0] -= TWO_PI * floorf(ex_phase[0] / TWO_PI + 0.5);

	for (m = 1; m <= model->L; m++)
	{
		/* generate excitation */
		if (model->voiced)
		{
			Ex[m].r = cosf(ex_phase[0] * m);
			Ex[m].i = sinf(ex_phase[0] * m);
		}
		else
		{
			/* When a few samples were tested I found that LPC filter
			   phase is not needed in the unvoiced case, but no harm in
			   keeping it.
			*/
			float phi = TWO_PI * (float)codec2_rand(&c2->next_rn) / CODEC2_RAND_MAX;
			Ex[m].r = cosf(phi);
			Ex[m].i = sinf(phi);
		}

		/* filter using LPC filter */
		A_[m].r = H[m].r * Ex[m].r - H[m].i * Ex[m].i;
		A_[m].i = H[m].i * Ex[m].r + H[m].r * Ex[m].i;

		/* modify sinusoidal phase */
		new_phi = atan2f(A_[m].i, A_[m].r + 1E-12);
		model->phi[m] = new_phi;
	}
}

void ear_protection(float *in_out, int n)
{
	float max_sample, over, gain;
	int i;

	/* find maximum sample in frame */
	max_sample = 0.0;
	for (i = 0; i < n; i++)
		if (in_out[i] > max_sample)
			max_sample = in_out[i];

	/* determine how far above set point */
	over = max_sample / 30000.0;

	/* If we are x dB over set point we reduce level by 2x dB, this
	   attenuates major excursions in amplitude (likely to be caused
	   by bit errors) more than smaller ones */
	if (over > 1.0)
	{
		gain = 1.0 / (over * over);
		for (i = 0; i < n; i++)
			in_out[i] *= gain;
	}
}

void analyse_one_frame(codec2_t *c2, model_t *model, const int16_t *speech)
{
	complex_t Sw[FFT_ENC];
	float pitch;
	int i;

	/* Read input speech */
	for (i = 0; i < M_PITCH - N_SAMP; i++)
		c2->Sn[i] = c2->Sn[i + N_SAMP];
	for (i = 0; i < N_SAMP; i++)
		c2->Sn[i + M_PITCH - N_SAMP] = speech[i];

	dft_speech(c2->fft_fwd_cfg, Sw, c2->Sn, c2->w);

	/* Estimate pitch */
	nlp(&c2->nlp, c2->Sn, N_SAMP, &pitch, &c2->prev_f0_enc);
	model->Wo = TWO_PI / pitch;
	model->L = M_PI / model->Wo;

	/* estimate model parameters */
	two_stage_pitch_refinement(model, Sw);

	/* estimate phases */
	estimate_amplitudes(model, Sw, 1);
	est_voicing_mbe(model, Sw, c2->W);
}

void synthesise_one_frame(codec2_t *c2, int16_t *speech, model_t *model,
						  complex_t *Aw, float gain)
{
	int i;

	/* LPC based phase synthesis */
	complex_t H[MAX_AMP + 1];
	sample_phase(model, H, Aw);
	phase_synth_zero_order(c2, model, &c2->ex_phase, H);
	postfilter(c2, model, &c2->bg_est);
	synthesise(c2->fftr_inv_cfg, c2->Sn_, model, c2->Pn, 1);

	for (i = 0; i < N_SAMP; i++)
	{
		c2->Sn_[i] *= gain;
	}

	ear_protection(c2->Sn_, N_SAMP);

	for (i = 0; i < N_SAMP; i++)
	{
		if (c2->Sn_[i] > 32767.0f)
			speech[i] = 32767;
		else if (c2->Sn_[i] < -32767.0f)
			speech[i] = -32767;
		else
			speech[i] = c2->Sn_[i];
	}
}

void pack(
	uint8_t *bitArray,		/* The output bit array. */
	unsigned int *bitIndex, /* Index into the string in BITS, not bytes.*/
	unsigned int field,		/* The bit field to be packed. */
	unsigned int fieldWidth /* Width of the field in BITS, not bytes. */
)
{
	/* Convert the field to Gray code, but only if the field is more than 1 bit */
	if (fieldWidth > 1)
		field = (field >> 1) ^ field;

	do
	{
		unsigned int bI = *bitIndex;
		unsigned int bitsLeft = 8 - (bI & 0x7);
		unsigned int sliceWidth = bitsLeft < fieldWidth ? bitsLeft : fieldWidth;
		unsigned int wordIndex = bI >> 3;

		bitArray[wordIndex] |= ((uint8_t)((field >> (fieldWidth - sliceWidth)) << (bitsLeft - sliceWidth)));

		*bitIndex = bI + sliceWidth;
		fieldWidth -= sliceWidth;
	} while (fieldWidth != 0);
}

int unpack(
	const uint8_t *bitArray, /* The input bit string. */
	unsigned int *bitIndex,	 /* Index into the string in BITS, not bytes.*/
	unsigned int fieldWidth	 /* Width of the field in BITS, not bytes. */
)
{
	unsigned int field = 0;
	unsigned int origWidth = fieldWidth;
	unsigned int t;

	do
	{
		unsigned int bI = *bitIndex;
		unsigned int bitsLeft = 8 - (bI & 0x7);
		unsigned int sliceWidth = bitsLeft < fieldWidth ? bitsLeft : fieldWidth;

		field |= (((bitArray[bI >> 3] >> (bitsLeft - sliceWidth)) &
				   ((1 << sliceWidth) - 1))
				  << (fieldWidth - sliceWidth));

		*bitIndex = bI + sliceWidth;
		fieldWidth -= sliceWidth;
	} while (fieldWidth != 0);

	if (origWidth > 1)
	{
		/* Convert from Gray code to binary. Works for maximum 8-bit fields. */
		t = field ^ (field >> 8);
		t ^= (t >> 4);
		t ^= (t >> 2);
		t ^= (t >> 1);
	}
	else
	{
		t = field;
	}

	return t;
}

int encode_Wo(float Wo, uint8_t bits)
{
	int index, Wo_levels = 1 << bits;
	float norm;

	norm = (Wo - W0_MIN) / (W0_MAX - W0_MIN);
	index = floorf(Wo_levels * norm + 0.5);
	if (index < 0)
		index = 0;
	if (index > (Wo_levels - 1))
		index = Wo_levels - 1;

	return index;
}

float decode_Wo(int index, int bits)
{
	float step;
	float Wo;
	int Wo_levels = 1 << bits;

	step = (W0_MAX - W0_MIN) / Wo_levels;
	Wo = W0_MIN + step * (index);

	return Wo;
}

void autocorrelate(float *Sn, /* frame of Nsam windowed speech samples */
				   float *Rn  /* array of P+1 autocorrelation coefficients */
)
{
	int i, j; /* loop variables */

	for (j = 0; j < LPC_ORD + 1; j++)
	{
		Rn[j] = 0.0;
		for (i = 0; i < M_PITCH - j; i++)
			Rn[j] += Sn[i] * Sn[i + j];
	}
}

void levinson_durbin(float *R,	 /* order+1 autocorrelation coeff */
					 float *lpcs /* order+1 LPC's */
)
{
	float a[LPC_ORD + 1][LPC_ORD + 1];
	float sum, e, k;
	int i, j; /* loop variables */

	e = R[0]; /* Equation 38a, Makhoul */

	for (i = 1; i <= LPC_ORD; i++)
	{
		sum = 0.0;
		for (j = 1; j <= i - 1; j++)
			sum += a[i - 1][j] * R[i - j];
		k = -1.0 * (R[i] + sum) / e; /* Equation 38b, Makhoul */
		if (fabsf(k) > 1.0)
			k = 0.0;

		a[i][i] = k;

		for (j = 1; j <= i - 1; j++)
			a[i][j] = a[i - 1][j] + k * a[i - 1][i - j]; /* Equation 38c, Makhoul */

		e *= (1 - k * k); /* Equation 38d, Makhoul */
	}

	for (i = 1; i <= LPC_ORD; i++)
		lpcs[i] = a[LPC_ORD][i];
	lpcs[0] = 1.0;
}

/*  float coef[]  	coefficients of the polynomial to be evaluated 	*/
/*  float x   		the point where polynomial is to be evaluated 	*/
/*  int order 		order of the polynomial 			            */
/*  NOTE: this function uses fully unrolled loop for LPC_ORD=10     */
inline float cheb_poly_eva(float *coef, float x)
{
	_Static_assert(LPC_ORD == 10, "cheb_poly_eva() assumes LPC_ORD=10");

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

/*  float *a 		     	lpc coefficients			*/
/*  float *freq 	      	LSP frequencies in radians  */
/*  NOTE: This function uses 6 bisections       		*/
int lpc_to_lsp(float *a, float *freq)
{
	float psuml, psumr, psumm, xl, xr, xm;
	int i, j, m;
	float *pt;	   /* ptr used for cheb_poly_eval(), whether P' or Q' */
	int roots = 0; /* number of roots found */
	float Q[LPC_ORD + 1];
	float P[LPC_ORD + 1];

	m = LPC_ORD / 2; /* order of P'(z) & Q'(z) polynimials 	*/

	P[0] = Q[0] = 1.0f;
	for (i = 1; i <= m; i++)
	{
		P[i] = a[i] + a[LPC_ORD + 1 - i] - P[i - 1];
		Q[i] = a[i] - a[LPC_ORD + 1 - i] + Q[i - 1];
	}

	for (i = 0; i < m; i++)
	{
		P[i] *= 2.0f;
		Q[i] *= 2.0f;
	}

	/* Search for a zero in P'(z) polynomial first and then alternate to Q'(z).
	Keep alternating between the two polynomials as each zero is found 	*/
	xl = 1.0; /* start at point xl = 1 		*/
	const float delta = LSP_DELTA1;

	for (j = 0; j < LPC_ORD; j++)
	{
		pt = (j & 1) ? Q : P; /* determines whether P' or Q' is eval. */

		xr = xl;
		psuml = cheb_poly_eva(pt, xl); /* evals poly. at xl 	*/

		while (xr >= -1.0f)
		{
			xr = xl - delta;			   /* interval spacing 	*/
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
	for (i = 0; i < LPC_ORD; i++)
	{
		freq[i] = acosf(freq[i]);
	}

	return (roots);
}

/*  float *freq         array of LSP frequencies in radians     	*/
/*  float *ak 		array of LPC coefficients 			*/
void lsp_to_lpc(float *lsp, float *ak)
{
	int i, j;
	float xout1, xout2, xin1, xin2;
	float *pw, *n1, *n2, *n3, *n4 = 0;
	float freq[LPC_ORD];
	float Wp[(LPC_ORD * 4) + 2];

	/* convert from radians to the x=cos(w) domain */
	for (i = 0; i < LPC_ORD; i++)
		freq[i] = cosf(lsp[i]);

	pw = Wp;

	/* initialise contents of array */
	for (i = 0; i <= 4 * (LPC_ORD / 2) + 1; i++)
	{ /* set contents of buffer to 0 */
		*pw++ = 0.0;
	}

	/* Set pointers up */
	pw = Wp;
	xin1 = 1.0;
	xin2 = 1.0;

	/* reconstruct P(z) and Q(z) by cascading second order polynomials
	  in form 1 - 2xz(-1) +z(-2), where x is the LSP coefficient */
	for (j = 0; j <= LPC_ORD; j++)
	{
		for (i = 0; i < (LPC_ORD / 2); i++)
		{
			n1 = pw + (i * 4);
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

void lpc_post_filter(float *Pw,
					 const float *A2,
					 const float *A2g,
					 int bass_boost,
					 float E)
{
	float e_before = 1E-4;
	float e_after = 1E-4;

	/* measure energy before post filtering */
	for (int i = 0; i < FFT_ENC / 2; i++)
		e_before += Pw[i];

	/* R(w) = |A_gamma| / |A| */
	for (int i = 0; i < FFT_ENC / 2; i++)
	{
		float R = sqrtf(A2g[i] / A2[i]);
		Pw[i] *= expf(LPCPF_TWO_BETA * logf(R));
		e_after += Pw[i];
	}

	/* normalise energy and apply LPC energy */
	float gain = e_before / e_after * E;

	for (int i = 0; i < FFT_ENC / 2; i++)
	{
		Pw[i] *= gain;
	}

	if (bass_boost)
	{
		/* add a slight boost to first 1 kHz to account for LP effect of PF */
		for (int i = 0; i < FFT_ENC / 8; i++)
		{
			Pw[i] *= 1.9327956f; // this value seems to maximize ViSQOL MOS
		}
	}
}

void aks_to_mag2(codec2_t *c2,
				 float *ak,		 /* LPC's */
				 model_t *model, /* sinusoidal model parameters for this frame */
				 float E,		 /* energy term */
				 float *snr,	 /* signal to noise ratio for this frame in dB */
				 int sim_pf,	 /* true to simulate a post filter */
				 int pf,		 /* true to enable actual LPC post filter */
				 int bass_boost, /* enable LPC filter 0-1kHz 3dB boost */
				 complex_t *Aw,	 /* output power spectrum */
				 float *A2)
{
	int i, m;	/* loop variables */
	int am, bm; /* limits of current band */
	float Em;	/* energy in band */
	float Am;	/* spectral amplitude sample */
	float signal, noise;

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
	complex_t Awg[FFT_ENC / 2 + 1];
	kiss_fftr(c2->fftr_fwd_cfg, ag, Awg);

	float A2g[FFT_ENC / 2];
	for (int i = 0; i < FFT_ENC / 2; i++)
	{
		A2g[i] = Awg[i].r * Awg[i].r + Awg[i].i * Awg[i].i + 1e-6f;
	}

	/* Determine power spectrum P(w) = E/(A(exp(jw))^2 ------------------------*/
	float Pw[FFT_ENC / 2];

	for (i = 0; i < FFT_ENC / 2; i++)
	{
		Pw[i] = 1.0f / A2[i];
	}

	if (pf)
		lpc_post_filter(Pw, A2, A2g, bass_boost, E);
	else
	{
		for (i = 0; i < FFT_ENC / 2; i++)
		{
			Pw[i] *= E;
		}
	}

	/* Determine magnitudes from P(w) ----------------------------------------*/
	/* when used just by decoder {A} might be all zeroes so init signal
	   and noise to prevent log(0) errors */
	signal = 1E-30;
	noise = 1E-32;

	for (m = 1; m <= model->L; m++)
	{
		am = (int)((m - 0.5) * model->Wo / FFT_R + 0.5);
		bm = (int)((m + 0.5) * model->Wo / FFT_R + 0.5);

		// FIXME: With arm_rfft_fast_f32 we have to use this
		// otherwise sometimes a to high bm is calculated
		// which causes trouble later in the calculation
		// chain
		// it seems for some reason model->Wo is calculated somewhat too high
		if (bm > FFT_ENC / 2)
		{
			bm = FFT_ENC / 2;
		}
		Em = 0.0;

		for (i = am; i < bm; i++)
			Em += Pw[i];
		Am = sqrtf(Em);

		signal += model->A[m] * model->A[m];
		noise += (model->A[m] - Am) * (model->A[m] - Am);

		/* This code significantly improves perf of LPC model, in
		   particular when combined with phase0.  The LPC spectrum tends
		   to track just under the peaks of the spectral envelope, and
		   just above nulls.  This algorithm does the reverse to
		   compensate - raising the amplitudes of spectral peaks, while
		   attenuating the null.  This enhances the formants, and
		   suppresses the energy between formants. */
		if (sim_pf)
		{
			if (Am > model->A[m])
				Am *= 0.7;
			if (Am < model->A[m])
				Am *= 1.4;
		}
		model->A[m] = Am;
	}

	*snr = 10.0 * log10f(signal / noise);
}

void apply_lpc_correction(model_t *model)
{
	if (model->Wo < (M_PI * 150.0f / 4000.0f))
	{
		model->A[1] *= 0.032f;
	}
}

float speech_to_uq_lsps(float *restrict lsp, float *restrict ak, float *restrict energy, const float *restrict Sn, const float *restrict w)
{
	int roots;
	float Wn[M_PITCH];
	float R[LPC_ORD + 1];
	float e, E;

	e = 0.0;
	for (int i = 0; i < M_PITCH; i++)
	{
		Wn[i] = Sn[i] * w[i];
		e += Wn[i] * Wn[i];
	}

	/* trap 0 energy case as LPC analysis will fail */
	if (e == 0.0f)
	{
		for (int i = 0; i < LPC_ORD; i++)
			lsp[i] = (M_PI / LPC_ORD) * (float)i;
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

int encode_energy(float e, int bits)
{
	int index, e_levels = 1 << bits;
	float e_min = E_MIN_DB;
	float e_max = E_MAX_DB;
	float norm;

	e = 10.0 * log10f(e);
	norm = (e - e_min) / (e_max - e_min);
	index = floorf(e_levels * norm + 0.5);
	if (index < 0)
		index = 0;
	if (index > (e_levels - 1))
		index = e_levels - 1;

	return index;
}

float decode_energy(int index, int bits)
{
	float step;
	float e;
	int e_levels = 1 << bits;

	step = (E_MAX_DB - E_MIN_DB) / e_levels;
	e = E_MIN_DB + step * (index);
	e = POW10F(e / 10.0f);

	return e;
}

/* float   *cb;	current VQ codebook		*/
/* float   *vec;	vector to quantise		*/
/* float   *w;         weighting vector                */
/* float   *se;		accumulated squared error 	*/
long quantise(const float *cb, float *vec, float *w, float *se)
{
	float e;	 /* current error		*/
	long besti;	 /* best index so far		*/
	float beste; /* best error so far		*/
	long j;
	int i;
	float diff;

	besti = 0;
	beste = 1E32;
	for (j = 0; j < 32; j++)
	{
		e = 0.0;
		for (i = 0; i < 1; i++)
		{
			diff = cb[j * 1 + i] - vec[i];
			e += diff * w[i] * diff * w[i]; // powf(diff*w[i],2.0);
		}
		if (e < beste)
		{
			beste = e;
			besti = j;
		}
	}

	*se += beste;

	return (besti);
}

void encode_lspds_scalar(uint16_t *indexes, float *lsp)
{
	int i;
	float lsp_hz[LPC_ORD];
	float lsp__hz[LPC_ORD];
	float dlsp[LPC_ORD];
	float dlsp_[LPC_ORD];
	float wt[LPC_ORD];
	const float *cb;
	float se;

	for (i = 0; i < LPC_ORD; i++)
	{
		wt[i] = 1.0;
	}

	/* convert from radians to Hz so we can use human readable
	   frequencies */
	for (i = 0; i < LPC_ORD; i++)
		lsp_hz[i] = (4000.0 / M_PI) * lsp[i];

	wt[0] = 1.0;
	for (i = 0; i < LPC_ORD; i++)
	{
		/* find difference from previous quantised lsp */
		if (i)
			dlsp[i] = lsp_hz[i] - lsp__hz[i - 1];
		else
			dlsp[0] = lsp_hz[0];

		cb = &delta_lsp_cb[i][0];
		indexes[i] = quantise(cb, &dlsp[i], wt, &se);
		dlsp_[i] = cb[indexes[i]];

		if (i)
			lsp__hz[i] = lsp__hz[i - 1] + dlsp_[i];
		else
			lsp__hz[0] = dlsp_[0];
	}
}

void decode_lspds_scalar(float *lsp_, int *indexes)
{
	int i;
	float lsp__hz[LPC_ORD];
	float dlsp_[LPC_ORD];
	const float *cb;

	for (i = 0; i < LPC_ORD; i++)
	{
		cb = &delta_lsp_cb[i][0];
		dlsp_[i] = cb[indexes[i]];

		if (i)
			lsp__hz[i] = lsp__hz[i - 1] + dlsp_[i];
		else
			lsp__hz[0] = dlsp_[0];

		lsp_[i] = (M_PI / 4000.0f) * lsp__hz[i];
	}
}

void interp_Wo(model_t *interp, /* interpolated model params */
			   model_t *prev,	/* previous frames model params                  */
			   model_t *next	/* next frames model params                      */
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

float interp_energy(float prev_e, float next_e)
{
	// return powf(10.0, (log10f(prev_e) + log10f(next_e))/2.0);
	return sqrtf(prev_e * next_e);
}

void interpolate_lsp(float *interp, float *prev, float *next)
{
	const float weight = 0.5f;
	int i;

	for (i = 0; i < LPC_ORD; i++)
		interp[i] = (1.0 - weight) * prev[i] + weight * next[i];
}

void codec2_encode(codec2_t *c2, uint8_t *bits, const int16_t *speech)
{
	model_t model;
	float ak[LPC_ORD + 1];
	float lsps[LPC_ORD];
	float e;
	int Wo_index, e_index;
	uint16_t lspd_indexes[LPC_ORD];
	int i;
	unsigned int nbit = 0;

	memset(bits, 0, 8); // 64 bits

	/* first 10ms analysis frame - we just want voicing */
	analyse_one_frame(c2, &model, speech);
	pack(bits, &nbit, model.voiced, 1);

	/* second 10ms analysis frame */
	analyse_one_frame(c2, &model, &speech[N_SAMP]);
	pack(bits, &nbit, model.voiced, 1);
	Wo_index = encode_Wo(model.Wo, WO_BITS);
	pack(bits, &nbit, Wo_index, WO_BITS);

	speech_to_uq_lsps(lsps, ak, &e, c2->Sn, c2->w);
	e_index = encode_energy(e, E_BITS);
	pack(bits, &nbit, e_index, E_BITS);

	encode_lspds_scalar(lspd_indexes, lsps);
	for (i = 0; i < LSPD_SCALAR_INDEXES; i++)
	{
		pack(bits, &nbit, lspd_indexes[i], 5); // codebook indices are 5 bits
	}
}

void codec2_decode(codec2_t *c2, int16_t *speech, const uint8_t *bits)
{
	model_t model[2];
	int lspd_indexes[LPC_ORD];
	float lsps[2][LPC_ORD];
	int Wo_index, e_index;
	float e[2];
	float snr;
	float ak[2][LPC_ORD + 1];
	int i, j;
	unsigned int nbit = 0;
	complex_t Aw[FFT_ENC];
	float A2[FFT_ENC / 2];

	/* only need to zero these out due to (unused) snr calculation */
	for (i = 0; i < 2; i++)
		for (j = 1; j <= MAX_AMP; j++)
			model[i].A[j] = 0.0;

	/* unpack bits from channel ------------------------------------*/

	/* this will partially fill the model params for the 2 x 10ms
	   frames */
	model[0].voiced = unpack(bits, &nbit, 1);
	model[1].voiced = unpack(bits, &nbit, 1);

	Wo_index = unpack(bits, &nbit, WO_BITS);
	model[1].Wo = decode_Wo(Wo_index, WO_BITS);
	model[1].L = M_PI / model[1].Wo;

	e_index = unpack(bits, &nbit, E_BITS);
	e[1] = decode_energy(e_index, E_BITS);

	for (i = 0; i < LSPD_SCALAR_INDEXES; i++)
	{
		lspd_indexes[i] = unpack(bits, &nbit, 5); // codebook indices are 5 bits
	}
	decode_lspds_scalar(&lsps[1][0], lspd_indexes);

	/* interpolate ------------------------------------------------*/
	/* Wo and energy are sampled every 20ms, so we interpolate just 1
	   10ms frame between 20ms samples */
	interp_Wo(&model[0], &c2->prev_model_dec, &model[1]);
	e[0] = interp_energy(c2->prev_e_dec, e[1]);

	/* LSPs are sampled every 20ms so we interpolate the frame in
	   between, then recover spectral amplitudes */
	interpolate_lsp(&lsps[0][0], c2->prev_lsps_dec, &lsps[1][0]);

	for (i = 0; i < 2; i++)
	{
		lsp_to_lpc(&lsps[i][0], &ak[i][0]);
		aks_to_mag2(c2, &ak[i][0], &model[i], e[i], &snr, 0,
					c2->lpc_pf, c2->bass_boost, Aw, A2);
		apply_lpc_correction(&model[i]);
		synthesise_one_frame(c2, &speech[N_SAMP * i], &model[i], Aw, 1.0f);
	}

	/* update memories for next frame ----------------------------*/
	c2->prev_model_dec = model[1];
	c2->prev_e_dec = e[1];
	for (i = 0; i < LPC_ORD; i++)
		c2->prev_lsps_dec[i] = lsps[1][i];
}

int main(void)
{
	struct timespec tick, tock;

	codec2_t c2;
	codec2_init(&c2);

	uint8_t encoded[8];
	int16_t speech[160];

	FILE *audio_in, *bitstream, *audio_out;

	audio_in = fopen("../../sample.raw", "rb");
	bitstream = fopen("../../modified_bitstream.bin", "wb");
	audio_out = fopen("../../modified_decoded.raw", "wb");

	clock_gettime(CLOCK_MONOTONIC, &tick);
	size_t rb;
	while ((rb = fread(speech, 1, sizeof(speech), audio_in)) > 0)
	{
		if (rb < sizeof(speech))
			memset(&speech[rb / 2], 0, sizeof(speech) - rb); //`break;` here instead to cut away the remaining audio (just like `c2enc` does)

		codec2_encode(&c2, encoded, speech);
		fwrite(encoded, sizeof(encoded), 1, bitstream);
	}
	clock_gettime(CLOCK_MONOTONIC, &tock);

	fclose(bitstream);

	long sec = tock.tv_sec - tick.tv_sec;
	long nsec = tock.tv_nsec - tick.tv_nsec;
	if (nsec < 0)
	{
		sec--;
		nsec += 1000000000L;
	}
	double elapsed_ms = sec * 1e3 + nsec * 1e-6;
	fprintf(stderr, "Elapsed time: %.3f ms\n", elapsed_ms);

	bitstream = fopen("../../modified_bitstream.bin", "rb");
	codec2_init(&c2);

	while (fread(encoded, 8, 1, bitstream) == 1)
	{
		codec2_decode(&c2, speech, encoded);
		fwrite(speech, sizeof(speech), 1, audio_out);
	}

	fclose(audio_in);
	fclose(bitstream);
	fclose(audio_out);

	/*for (uint8_t i = 0; i < 160; i++)
		speech[i] = 0.5f * sinf(i / 80.0f * TWO_PI);

	for (uint8_t j = 0; j < 10; j++)
	{
		codec2_encode(&c2, encoded, speech);

		for (uint8_t i = 0; i < 8; i++)
			fprintf(stderr, "%02X ", encoded[i]);
		fprintf(stderr, "\n");
	}*/

	return 0;
}
