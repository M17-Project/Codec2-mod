/*
  This work uses code written by David Rowe VK5DGR et al.
  https://github.com/drowe67/codec2
*/

#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "kiss_fft.h"
#include "kiss_fftr.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N_S 0.01   /* internal proc frame length in secs   */
#define TW_S 0.005 /* trapezoidal synth window overlap     */
#define MAX_AMP 80 /* maximum number of harmonics          */

#define FFT_ENC 512	 /* size of FFT used for encoder         */
#define FFT_DEC 512	 /* size of FFT used in decoder          */
#define V_THRESH 6.0 /* voicing threshold in dB              */
#define LPC_ORD 10	 /* LPC order                            */

/* Pitch estimation defines */
#define M_PITCH_S 0.0400f /* pitch analysis window in s           */
#define P_MIN_S 0.0025f	  /* minimum pitch period in s            */
#define P_MAX_S 0.0200f	  /* maximum pitch period in s            */

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
#define T 0.1			/* threshold for local minima candidate */
#define F0_MAX 500
#define CNLP 0.3	/* post processor constant              */
#define NLP_NTAP 48 /* Decimation LPF order */

/* quantizers & LPC */
#define WO_BITS 7
#define WO_LEVELS (1 << WO_BITS)
#define WO_DT_BITS 3

#define E_BITS 5
#define E_LEVELS (1 << E_BITS)
#define E_MIN_DB -10.0
#define E_MAX_DB 40.0

#define LSP_SCALAR_INDEXES 10
#define LSPD_SCALAR_INDEXES 10
#define LSP_PRED_VQ_INDEXES 3

#define WO_E_BITS 8

#define LPCPF_GAMMA 0.5
#define LPCPF_BETA 0.2

#define LSP_DELTA1 0.01 /* grid spacing for LSP root searches */

/* bandpass filter taps? */
#define BPF_N 101

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
	int m;
	float w[PMAX_M / DEC];	 /* DFT window                   */
	float sq[PMAX_M];		 /* squared speech samples       */
	float mem_x, mem_y;		 /* memory for notch filter      */
	float mem_fir[NLP_NTAP]; /* decimation FIR filter memory */
	kiss_fft_cfg fft_cfg;	 /* kiss FFT config              */
	float *Sn16k;			 /* Fs=16kHz input speech vector */
} nlp_t;

typedef struct codec2_t
{
	kiss_fft_cfg fft_fwd_cfg;		   /* forward FFT config                        */
	kiss_fftr_cfg fftr_fwd_cfg;		   /* forward real FFT config                   */
	float w[M_PITCH];				   /* [m_pitch] time domain hamming window      */
	float W[FFT_ENC];				   /* DFT of w[]                                */
	float Pn[2 * N_SAMP];			   /* [2*n_samp] trapezoidal synthesis window   */
	float bpf_buf[BPF_N + 4 * N_SAMP]; /* buffer for band pass filter               */
	float Sn[M_PITCH];				   /* [m_pitch] input speech                    */
	float hpf_states[2];			   /* high pass filter states                   */
	nlp_t nlp;						   /* pitch predictor states                    */

	kiss_fftr_cfg fftr_inv_cfg;	  /* inverse FFT config                        */
	float Sn_[2 * N_SAMP];		  /* [2*n_samp] synthesised output speech      */
	float ex_phase;				  /* excitation model phase track              */
	float bg_est;				  /* background noise estimate for post filter */
	float prev_f0_enc;			  /* previous frame's f0    estimate           */
	model_t prev_model_dec;		  /* previous frame's model parameters         */
	float prev_lsps_dec[LPC_ORD]; /* previous frame's LSPs                     */
	float prev_e_dec;			  /* previous frame's LPC energy               */

	int lpc_pf;		/* LPC post filter on                        */
	int bass_boost; /* LPC post filter bass boost                */
	float beta;		/* LPC post filter parameters                */
	float gamma;

	float xq_enc[2]; /* joint pitch and energy VQ states          */
	float xq_dec[2];

	/* newamp1 states */
	// float rate_K_sample_freqs_kHz[NEWAMP1_K];
	// float prev_rate_K_vec_[NEWAMP1_K];
	float Wo_left;
	int voicing_left;
	kiss_fft_cfg phase_fft_fwd_cfg;
	kiss_fft_cfg phase_fft_inv_cfg;
	float se; /* running sum of squared error */
	int nse;  /* number of terms in sum       */

	bool post_filter_en;
	bool eq_en;
} codec2_t;

void make_analysis_window(kiss_fft_cfg fft_fwd_cfg, float *w, float *W)
{
	float m;
	complex_t wshift[FFT_ENC];
	int i, j;

	m = 0.0;
	for (i = 0; i < M_PITCH / 2 - NW / 2; i++)
		w[i] = 0.0f;
	for (i = M_PITCH / 2 - NW / 2, j = 0; i < M_PITCH / 2 + NW / 2; i++, j++)
	{
		w[i] = 0.5f - 0.5f * cosf(2.0f * M_PI * j / (NW - 1));
		m += w[i] * w[i];
	}
	for (i = M_PITCH / 2 + NW / 2; i < M_PITCH; i++)
		w[i] = 0.0f;

	m = 1.0f / sqrtf(m * FFT_ENC);
	for (i = 0; i < M_PITCH; i++)
	{
		w[i] *= m;
	}

	complex_t temp[FFT_ENC];

	for (i = 0; i < FFT_ENC; i++)
	{
		wshift[i].r = 0.0;
		wshift[i].i = 0.0;
	}
	for (i = 0; i < NW / 2; i++)
		wshift[i].r = w[i + M_PITCH / 2];
	for (i = FFT_ENC - NW / 2, j = M_PITCH / 2 - NW / 2; i < FFT_ENC; i++, j++)
		wshift[i].r = w[j];

	kiss_fft(fft_fwd_cfg, wshift, temp);

	for (i = 0; i < FFT_ENC / 2; i++)
	{
		W[i] = temp[i + FFT_ENC / 2].r;
		W[i + FFT_ENC / 2] = temp[i].r;
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
	int i;
	int m = M_PITCH;

	nlp->m = m;

	assert(m <= PMAX_M);

	for (i = 0; i < m / DEC; i++)
	{
		nlp->w[i] = 0.5 - 0.5 * cosf(2.0f * M_PI * i / (m / DEC - 1));
	}

	for (i = 0; i < PMAX_M; i++)
		nlp->sq[i] = 0.0;

	nlp->mem_x = 0.0;
	nlp->mem_y = 0.0;

	for (i = 0; i < NLP_NTAP; i++)
		nlp->mem_fir[i] = 0.0;

	nlp->fft_cfg = kiss_fft_alloc(PE_FFT_SIZE, 0, NULL, NULL);
	assert(nlp->fft_cfg != NULL);
}

void codec2_init(codec2_t *c2)
{
	int i, l;

	for (i = 0; i < M_PITCH; i++)
		c2->Sn[i] = 1.0f;

	c2->hpf_states[0] = c2->hpf_states[1] = 0.0;

	for (i = 0; i < 2 * N_SAMP; i++)
		c2->Sn_[i] = 0.0f;

	c2->fft_fwd_cfg = kiss_fft_alloc(FFT_ENC, 0, NULL, NULL);
	c2->fftr_fwd_cfg = kiss_fftr_alloc(FFT_ENC, 0, NULL, NULL);

	make_analysis_window(c2->fft_fwd_cfg, c2->w, c2->W);
	make_synthesis_window(c2->Pn);

	c2->fftr_inv_cfg = kiss_fftr_alloc(FFT_DEC, 1, NULL, NULL);
	c2->prev_f0_enc = 1.0f / P_MAX_S;
	c2->bg_est = 0.0;
	c2->ex_phase = 0.0;

	for (l = 1; l <= MAX_AMP; l++)
		c2->prev_model_dec.A[l] = 0.0;
	c2->prev_model_dec.Wo = 2.0f * M_PI / P_MAX;
	c2->prev_model_dec.L = M_PI / c2->prev_model_dec.Wo;
	c2->prev_model_dec.voiced = 0;

	for (i = 0; i < LPC_ORD; i++)
	{
		c2->prev_lsps_dec[i] = i * M_PI / (LPC_ORD + 1);
	}
	c2->prev_e_dec = 1;

	nlp_init(&c2->nlp);

	c2->lpc_pf = 1;
	c2->bass_boost = 1;
	c2->beta = LPCPF_BETA;
	c2->gamma = LPCPF_GAMMA;

	c2->xq_enc[0] = c2->xq_enc[1] = 0.0;
	c2->xq_dec[0] = c2->xq_dec[1] = 0.0;

	c2->se = 0.0;
	c2->nse = 0;
	c2->post_filter_en = true;

	for (i = 0; i < BPF_N + 4 * N_SAMP; i++)
		c2->bpf_buf[i] = 0.0;
}

void codec2_fft_inplace(kiss_fft_cfg cfg, kiss_fft_cpx *inout)
{
	kiss_fft_cpx in[FFT_ENC];
	memcpy(in, inout, FFT_ENC * sizeof(kiss_fft_cpx));
	kiss_fft(cfg, in, inout);
}

void dft_speech(kiss_fft_cfg fft_fwd_cfg, complex_t *Sw, float *Sn, float *w)
{
	int i;

	for (i = 0; i < FFT_ENC; i++)
	{
		Sw[i].r = 0.0f;
		Sw[i].i = 0.0f;
	}

	/* Centre analysis window on time axis, we need to arrange input
	   to FFT this way to make FFT phases correct */
	/* move 2nd half to start of FFT input vector */
	for (i = 0; i < NW / 2; i++)
		Sw[i].r = Sn[i + M_PITCH / 2] * w[i + M_PITCH / 2];

	/* move 1st half to end of FFT input vector */
	for (i = 0; i < NW / 2; i++)
		Sw[FFT_ENC - NW / 2 + i].r =
			Sn[i + M_PITCH / 2 - NW / 2] * w[i + M_PITCH / 2 - NW / 2];

	codec2_fft_inplace(fft_fwd_cfg, Sw);
}

void estimate_amplitudes(model_t *model, complex_t *Sw, int est_phase)
{
	int i, m;	/* loop variables */
	int am, bm; /* bounds of current harmonic */
	float den;	/* denominator of amplitude expression */

	float r = 2.0f * M_PI / FFT_ENC;
	float one_on_r = 1.0f / r;

	for (m = 1; m <= model->L; m++)
	{
		/* Estimate ampltude of harmonic */
		den = 0.0f;
		am = (int)((m - 0.5f) * model->Wo * one_on_r + 0.5f);
		bm = (int)((m + 0.5f) * model->Wo * one_on_r + 0.5f);

		for (i = am; i < bm; i++)
		{
			den += Sw[i].r * Sw[i].r + Sw[i].i * Sw[i].i;
		}

		model->A[m] = sqrtf(den);

		if (est_phase)
		{
			int b = (int)(m * model->Wo / r + 0.5); /* DFT bin of centre of current harmonic */

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
	int m;			   /* loop variable */
	int b;			   /* bin for current harmonic centre */
	float E;		   /* energy for current pitch*/
	float Wo;		   /* current "test" fundamental freq. */
	float Wom;		   /* Wo that maximises E */
	float Em;		   /* mamimum energy */
	float r, one_on_r; /* number of rads/bin */
	float p;		   /* current pitch */

	/* Initialisation */
	model->L = M_PI / model->Wo; /* use initial pitch est. for L */
	Wom = model->Wo;
	Em = 0.0;
	r = 2.0f * M_PI / FFT_ENC;
	one_on_r = 1.0 / r;

	/* Determine harmonic sum for a range of Wo values */
	for (p = pmin; p <= pmax; p += pstep)
	{
		E = 0.0;
		Wo = 2.0f * M_PI / p;

		float bFloat = Wo * one_on_r;
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
	pmax = 2.0f * M_PI / model->Wo + 5;
	pmin = 2.0f * M_PI / model->Wo - 5;
	pstep = 1.0;
	hs_pitch_refinement(model, Sw, pmin, pmax, pstep);

	/* Fine refinement */
	pmax = 2.0f * M_PI / model->Wo + 1;
	pmin = 2.0f * M_PI / model->Wo - 1;
	pstep = 0.25;
	hs_pitch_refinement(model, Sw, pmin, pmax, pstep);

	/* Limit range */
	if (model->Wo < 2.0f * M_PI / P_MAX)
		model->Wo = 2.0f * M_PI / P_MAX;
	if (model->Wo > 2.0f * M_PI / P_MIN)
		model->Wo = 2.0f * M_PI / P_MIN;

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
		al = ceilf((l - 0.5) * Wo * FFT_ENC / (2.0f * M_PI));
		bl = ceilf((l + 0.5) * Wo * FFT_ENC / (2.0f * M_PI));

		/* Estimate amplitude of harmonic assuming harmonic is totally voiced */

		offset = FFT_ENC / 2 - l * Wo * FFT_ENC / (2.0f * M_PI) + 0.5;
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
		sixty = 60.0f * 2.0f * M_PI / SAMP_RATE;
		if ((eratio < -4.0f) && (model->Wo <= sixty))
			model->voiced = 0;
	}

	return snr;
}

float nlp(
	void *nlp_state, float *Sn, /* input speech vector */
	int n,						/* frames shift (no. new samples in Sn[])             */
	float *pitch,				/* estimated pitch period in samples at current Fs    */
	float *prev_f0				/* previous pitch f0 in Hz, memory for pitch tracking */
)
{
	nlp_t *nlp;
	float notch;			   /* current notch filter output          */
	complex_t Fw[PE_FFT_SIZE]; /* DFT of squared signal (input/output) */
	float gmax;
	int gmax_bin;
	int m, i, j;
	float best_f0;

	assert(nlp_state != NULL);
	nlp = (nlp_t *)nlp_state;
	m = nlp->m;

	/* Square, notch filter at DC, and LP filter vector */

	/* Square latest input samples */
	for (i = m - n; i < m; i++)
	{
		nlp->sq[i] = Sn[i] * Sn[i];
	}

	for (i = m - n; i < m; i++)
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
	for (i = m - n; i < m; i++)
	{

		for (j = 0; j < NLP_NTAP - 1; j++)
			nlp->mem_fir[j] = nlp->mem_fir[j + 1];
		nlp->mem_fir[NLP_NTAP - 1] = nlp->sq[i];

		nlp->sq[i] = 0.0;
		for (j = 0; j < NLP_NTAP; j++)
			nlp->sq[i] += nlp->mem_fir[j] * nlp_fir[j];
	}

	/* Decimate and DFT */
	for (i = 0; i < PE_FFT_SIZE; i++)
	{
		Fw[i].r = 0.0;
		Fw[i].i = 0.0;
	}
	for (i = 0; i < m / DEC; i++)
	{
		Fw[i].r = nlp->sq[i * DEC] * nlp->w[i];
	}

	// FIXME: check if this can be converted to a real fft
	// since all imag inputs are 0
	codec2_fft_inplace(nlp->fft_cfg, Fw);

	for (i = 0; i < PE_FFT_SIZE; i++)
		Fw[i].r = Fw[i].r * Fw[i].r + Fw[i].i * Fw[i].i;

	/* todo: express everything in f0, as pitch in samples is dep on Fs */
	int pmin = P_MIN;
	int pmax = P_MAX;

	/* find global peak */
	gmax = 0.0f;
	gmax_bin = PE_FFT_SIZE * DEC / pmax;
	for (i = PE_FFT_SIZE * DEC / pmax; i <= PE_FFT_SIZE * DEC / pmin; i++)
	{
		if (Fw[i].r > gmax)
		{
			gmax = Fw[i].r;
			gmax_bin = i;
		}
	}

	best_f0 = post_process_sub_multiples(Fw, gmax, gmax_bin, prev_f0);

	/* Shift samples in buffer to make room for new samples */
	for (i = 0; i < m - n; i++)
		nlp->sq[i] = nlp->sq[i + n];

	/* return pitch period in samples and F0 estimate */
	*pitch = (float)SAMP_RATE / best_f0;

	*prev_f0 = best_f0;

	return (best_f0);
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
	model->Wo = 2.0f * M_PI / pitch;
	model->L = M_PI / model->Wo;

	/* estimate model parameters */
	two_stage_pitch_refinement(model, Sw);

	/* estimate phases when doing ML experiments */
	estimate_amplitudes(model, Sw, 0);
	est_voicing_mbe(model, Sw, c2->W);
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

void autocorrelate(float Sn[], /* frame of Nsam windowed speech samples */
				   float Rn[]  /* array of P+1 autocorrelation coefficients */
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

void levinson_durbin(float R[],	  /* order+1 autocorrelation coeff */
					 float lpcs[] /* order+1 LPC's */
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

static float cheb_poly_eva(float *coef, float x)
/*  float coef[]  	coefficients of the polynomial to be evaluated 	*/
/*  float x   		the point where polynomial is to be evaluated 	*/
/*  int order 		order of the polynomial 			*/
{
	int i;
	float *t, *u, *v, sum;
	float Tm[(LPC_ORD / 2) + 1];

	/* Initialise pointers */

	t = Tm; /* T[i-2] 			*/
	*t++ = 1.0;
	u = t--; /* T[i-1] 			*/
	*u++ = x;
	v = u--; /* T[i] 			*/

	/* Evaluate chebyshev series formulation using iterative approach 	*/

	for (i = 2; i <= LPC_ORD / 2; i++)
		*v++ = (2 * x) * (*u++) - *t++; /* T[i] = 2*x*T[i-1] - T[i-2]	*/

	sum = 0.0; /* initialise sum to zero 	*/
	t = Tm;	   /* reset pointer 		*/

	/* Evaluate polynomial and return value also free memory space */

	for (i = 0; i <= LPC_ORD / 2; i++)
		sum += coef[(LPC_ORD / 2) - i] * *t++;

	return sum;
}

/*  float *a 		     	lpc coefficients			*/
/*  float *freq 	      	LSP frequencies in radians      	*/
/*  int nb			number of sub-intervals (4) 		*/
int lpc_to_lsp(float *a, float *freq, int nb)
{
	float psuml, psumr, psumm, temp_xr, xl, xr, xm = 0;
	float temp_psumr;
	int i, j, m, flag, k;
	float *px; /* ptrs of respective P'(z) & Q'(z)	*/
	float *qx;
	float *p;
	float *q;
	float *pt;	   /* ptr used for cheb_poly_eval()
					  whether P' or Q' 			*/
	int roots = 0; /* number of roots found 	        */
	float Q[LPC_ORD + 1];
	float P[LPC_ORD + 1];

	flag = 1;
	m = LPC_ORD / 2; /* order of P'(z) & Q'(z) polynimials 	*/

	/* Allocate memory space for polynomials */
	/* determine P'(z)'s and Q'(z)'s coefficients where
	  P'(z) = P(z)/(1 + z^(-1)) and Q'(z) = Q(z)/(1-z^(-1)) */
	px = P; /* initilaise ptrs */
	qx = Q;
	p = px;
	q = qx;
	*px++ = 1.0;
	*qx++ = 1.0;
	for (i = 1; i <= m; i++)
	{
		*px++ = a[i] + a[LPC_ORD + 1 - i] - *p++;
		*qx++ = a[i] - a[LPC_ORD + 1 - i] + *q++;
	}
	px = P;
	qx = Q;
	for (i = 0; i < m; i++)
	{
		*px = 2 * *px;
		*qx = 2 * *qx;
		px++;
		qx++;
	}
	px = P; /* re-initialise ptrs 			*/
	qx = Q;

	/* Search for a zero in P'(z) polynomial first and then alternate to Q'(z).
	Keep alternating between the two polynomials as each zero is found 	*/
	xr = 0;	  /* initialise xr to zero 		*/
	xl = 1.0; /* start at point xl = 1 		*/

	for (j = 0; j < LPC_ORD; j++)
	{
		if (j % 2) /* determines whether P' or Q' is eval. */
			pt = qx;
		else
			pt = px;

		psuml = cheb_poly_eva(pt, xl); /* evals poly. at xl 	*/
		flag = 1;
		while (flag && (xr >= -1.0))
		{
			xr = xl - LSP_DELTA1;		   /* interval spacing 	*/
			psumr = cheb_poly_eva(pt, xr); /* poly(xl-delta_x) 	*/
			temp_psumr = psumr;
			temp_xr = xr;

			/* if no sign change increment xr and re-evaluate
			   poly(xr). Repeat til sign change.  if a sign change has
			   occurred the interval is bisected and then checked again
			   for a sign change which determines in which interval the
			   zero lies in.  If there is no sign change between poly(xm)
			   and poly(xl) set interval between xm and xr else set
			   interval between xl and xr and repeat till root is located
			   within the specified limits  */
			if (((psumr * psuml) < 0.0) || (psumr == 0.0))
			{
				roots++;

				psumm = psuml;
				for (k = 0; k <= nb; k++)
				{
					xm = (xl + xr) / 2; /* bisect the interval 	*/
					psumm = cheb_poly_eva(pt, xm);
					if (psumm * psuml > 0.)
					{
						psuml = psumm;
						xl = xm;
					}
					else
					{
						psumr = psumm;
						xr = xm;
					}
				}

				/* once zero is found, reset initial interval to xr 	*/
				freq[j] = (xm);
				xl = xm;
				flag = 0; /* reset flag for next search 	*/
			}
			else
			{
				psuml = temp_psumr;
				xl = temp_xr;
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

float speech_to_uq_lsps(float *lsp, float *ak, float *Sn, float *w)
{
	int i, roots;
	float Wn[M_PITCH];
	float R[LPC_ORD + 1];
	float e, E;

	e = 0.0;
	for (i = 0; i < M_PITCH; i++)
	{
		Wn[i] = Sn[i] * w[i];
		e += Wn[i] * Wn[i];
	}

	/* trap 0 energy case as LPC analysis will fail */

	if (e == 0.0)
	{
		for (i = 0; i < LPC_ORD; i++)
			lsp[i] = (M_PI / LPC_ORD) * (float)i;
		return 0.0;
	}

	autocorrelate(Wn, R);
	levinson_durbin(R, ak);

	E = 0.0;
	for (i = 0; i <= LPC_ORD; i++)
		E += ak[i] * R[i];

	/* 15 Hz BW expansion as I can't hear the difference and it may help
	   help occasional fails in the LSP root finding.  Important to do this
	   after energy calculation to avoid -ve energy values.
	*/

	for (i = 0; i <= LPC_ORD; i++)
		ak[i] *= powf(0.994, (float)i);

	roots = lpc_to_lsp(ak, lsp, 5);
	if (roots != LPC_ORD)
	{
		/* if root finding fails use some benign LSP values instead */
		for (i = 0; i < LPC_ORD; i++)
			lsp[i] = (M_PI / LPC_ORD) * (float)i;
	}

	return E;
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

	e = speech_to_uq_lsps(lsps, ak, c2->Sn, c2->w);
	e_index = encode_energy(e, E_BITS);
	pack(bits, &nbit, e_index, E_BITS);

	encode_lspds_scalar(lspd_indexes, lsps);
	for (i = 0; i < LSPD_SCALAR_INDEXES; i++)
	{
		pack(bits, &nbit, lspd_indexes[i], 5); // every codebook entry is 5 bits long
	}
}

int main(void)
{
	codec2_t c2;
	codec2_init(&c2);

	uint8_t encoded[8] = {0};
	int16_t speech[160];

	for (uint8_t i = 0; i < 160; i++)
		speech[i] = 0.5f * sinf(i / 80.0f * 2.0f * M_PI);

	for (uint8_t j = 0; j < 10; j++)
	{
		codec2_encode(&c2, encoded, speech);

		for (uint8_t i = 0; i < 8; i++)
			fprintf(stderr, "%02X ", encoded[i]);
		fprintf(stderr, "\n");
	}

	return 0;
}
