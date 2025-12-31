#include "codec2_internal.h"
#include "analysis.h"
#include "nlp.h"
#include "util.h"

static void make_analysis_window(kiss_fft_cfg fft_fwd_cfg, float *w, float *W)
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

static void hs_pitch_refinement(model_t *restrict model, const complex_t *restrict Sw, float pmin, float pmax, float pstep)
{
	int b;	   /* bin for current harmonic centre */
	float E;   /* energy for current pitch*/
	float Wo;  /* current "test" fundamental freq. */
	float Wom; /* Wo that maximises E */
	float Em;  /* mamimum energy */

	/* Initialisation */
	model->L = M_PI / model->Wo; /* use initial pitch est. for L */
	Wom = model->Wo;
	Em = 0.0;

	/* Determine harmonic sum for a range of Wo values */
	for (float p = pmin; p <= pmax; p += pstep)
	{
		E = 0.0;
		Wo = TWO_PI / p;

		float bFloat = Wo * FFT_1_R;
		float currentBFloat = bFloat;

		/* Sum harmonic magnitudes */
		for (int m = 1; m <= model->L; m++)
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

static void two_stage_pitch_refinement(model_t *restrict model, const complex_t *restrict Sw)
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

static void estimate_amplitudes(model_t *model, const complex_t *Sw, int est_phase)
{
	for (int m = 1; m <= model->L; m++)
	{
		/* Estimate ampltude of harmonic */
		float den = 0.0f; /* denominator of amplitude expression */

		/* bounds of current harmonic */
		int am = (int)((m - 0.5f) * model->Wo * FFT_1_R + 0.5f);
		int bm = (int)((m + 0.5f) * model->Wo * FFT_1_R + 0.5f);

		// clamp
		if (am < 0)
			am = 0;
		if (bm > FFT_ENC / 2)
			bm = FFT_ENC / 2;

		for (int i = am; i < bm; i++)
		{
			den += Sw[i].r * Sw[i].r + Sw[i].i * Sw[i].i;
		}

		model->A[m] = sqrtf(den);

		/* recompute phases only for voiced speech :-) */
		if (est_phase && model->voiced)
		{
			int b = (int)(m * model->Wo / FFT_R + 0.5); /* DFT bin of centre of current harmonic */
			if (b >= FFT_ENC / 2)
				b = FFT_ENC / 2 - 1;

			/* Estimate phase of harmonic, this is expensive in CPU for
			   embedded devices, so we make it an option */
			model->phi[m] = fast_atan2f(Sw[b].i, Sw[b].r);
		}
	}
}

static float est_voicing_mbe(model_t *restrict model, const complex_t *restrict Sw, const float *restrict W)
{
	complex_t Am; /* amplitude sample for this band */
	int offset;	  /* centers Hw[] about current harmonic */
	float den;	  /* denominator of Am expression */
	float error;  /* accumulated error between original and synthesised */
	float Wo;
	float sig, snr;
	float elow, ehigh, eratio;
	float sixty;
	complex_t Ew;

	Ew.r = 0;
	Ew.i = 0;

	int l_1000hz = model->L * 1000.0 / (SAMP_RATE / 2);
	sig = 1E-4;
	for (int l = 1; l <= l_1000hz; l++)
	{
		sig += model->A[l] * model->A[l];
	}

	Wo = model->Wo;
	error = 1E-4;

	/* Just test across the harmonics in the first 1000 Hz */
	for (int l = 1; l <= l_1000hz; l++)
	{
		Am.r = 0.0;
		Am.i = 0.0;
		den = 0.0;

		int al = ceilf((l - 0.5) * Wo * FFT_ENC / TWO_PI);
		int bl = ceilf((l + 0.5) * Wo * FFT_ENC / TWO_PI);

		/* Estimate amplitude of harmonic assuming harmonic is totally voiced */
		offset = FFT_ENC / 2 - l * Wo * FFT_ENC / TWO_PI + 0.5;

		for (int m = al; m < bl; m++)
		{
			int idx = offset + m;
			if ((unsigned)idx < FFT_ENC)
			{
				Am.r += Sw[m].r * W[idx];
				Am.i += Sw[m].i * W[idx];
				den += W[idx] * W[idx];
			}
		}

		Am.r = Am.r / den;
		Am.i = Am.i / den;

		/* Determine error between estimated harmonic and original */
		for (int m = al; m < bl; m++)
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
	for (int l = 1; l <= l_2000hz; l++)
	{
		elow += model->A[l] * model->A[l];
	}
	for (int l = l_2000hz; l <= l_4000hz; l++)
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

static void dft_speech(kiss_fft_cfg fft_fwd_cfg, complex_t *Sw, const float *Sn, const float *w)
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

void analyse_one_frame(
	codec2_t *c2,
	model_t *model,
	const int16_t *speech)
{
	complex_t *Sw = c2->fft_buffer; // reuse scratch array
	float pitch;

	/* Read input speech */
	for (int i = 0; i < M_PITCH - N_SAMP; i++)
		c2->Sn[i] = c2->Sn[i + N_SAMP];
	for (int i = 0; i < N_SAMP; i++)
		c2->Sn[i + M_PITCH - N_SAMP] = speech[i];

	dft_speech(c2->fft_fwd_cfg, Sw, c2->Sn, c2->w);

	/* Estimate pitch */
	nlp(&c2->nlp, c2->Sn, &pitch, &c2->prev_f0_enc);
	model->Wo = TWO_PI / pitch;
	model->L = M_PI / model->Wo;

	/* estimate model parameters */
	two_stage_pitch_refinement(model, Sw);

	/* estimate phases */
	estimate_amplitudes(model, Sw, 1);
	est_voicing_mbe(model, Sw, c2->W);
}

void analysis_init(codec2_t *c2)
{
	make_analysis_window(c2->fft_fwd_cfg, c2->w, c2->W);
}
