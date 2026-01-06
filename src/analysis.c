#include "codec2_internal.h"
#include "analysis.h"
#include "nlp.h"
#include "util.h"

static void make_analysis_window(
	codec2_t *c2,
	kiss_fft_cfg fft_fwd_cfg,
	float *restrict w,
	float *restrict W)
{
	static const int mp2 = M_PITCH / 2;
	static const int nw2 = NW / 2;
	static const int fe2 = FFT_ENC / 2;

	complex_t *restrict wshift = c2->fft_buffer;

	/* zero time-domain window */
	memset(w, 0, M_PITCH * sizeof(float));

	/* build Hann window */
	float m = 0.0f;

	for (int j = 0; j < NW; j++)
	{
		float phase = TWO_PI * j / (NW - 1);
		float val = 0.5f - 0.5f * cosf(phase);
		int i = mp2 - nw2 + j;
		w[i] = val;
		m += val * val;
	}

	/* normalize (bit-exact with original) */
	const float scale = 1.0f / sqrtf(m * FFT_ENC);
	for (int i = 0; i < M_PITCH; i++)
		w[i] *= scale;

	/* zero FFT buffer */
	memset(wshift, 0, FFT_ENC * sizeof(*wshift));

	/* circular shift */
	for (int i = 0; i < nw2; i++)
		wshift[i].r = w[i + mp2];

	for (int i = FFT_ENC - nw2, j = mp2 - nw2; i < FFT_ENC; i++, j++)
		wshift[i].r = w[j];

	/* FFT */
	kiss_fft(fft_fwd_cfg, wshift, wshift);

	/* rearrange frequency response */
	for (int i = 0; i < fe2; i++)
	{
		W[i] = wshift[i + fe2].r;
		W[i + fe2] = wshift[i].r;
	}
}

static void hs_pitch_refinement(model_t *restrict model, const complex_t *restrict Sw, float pmin, float pmax, float pstep)
{
	float E;   /* energy for current pitch*/
	float Wo;  /* current "test" fundamental freq. */
	float Wom; /* Wo that maximises E */
	float Em;  /* mamimum energy */

	/* Initialisation */
	model->L = M_PI / model->Wo; /* use initial pitch est. for L */
	Wom = model->Wo;
	Em = 0.0f;

	/* Determine harmonic sum for a range of Wo values */
	for (float p = pmin; p <= pmax; p += pstep)
	{
		E = 0.0f;
		Wo = TWO_PI / p;

		float bFloat = Wo * FFT_1_R;
		float currentB = bFloat; // bin for current harmonic centre

		// limiting m to L=M_PI/Wo causes problems
		// we test the same number of harmonics per candidate
		for (int m = 1; m <= model->L; m++)
		{
			int b = (int)(currentB + 0.5f);

			// early breaking causes problems
			// reading past Nyquist adds some noise that helps
			// if (b >= FFT_DEC / 2)
			//	break;

			E += Sw[b].r * Sw[b].r + Sw[b].i * Sw[b].i;
			currentB += bFloat;
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

static void est_voicing_mbe(model_t *restrict model, const complex_t *restrict Sw, const float *restrict W)
{
	const float Wo = model->Wo;
	const float Wo_bin = Wo * FFT_ENC / TWO_PI;

	/*
	   Determine the ratio of low frequency to high frequency energy,
	   voiced speech tends to be dominated by low frequency energy,
	   unvoiced by high frequency. This measure can be used to
	   determine if we have made any gross errors.
	*/
	const int l_1000hz = model->L * 1000.0f / (SAMP_RATE / 2);
	const int l_2000hz = model->L * 2000.0f / (SAMP_RATE / 2);
	const int l_4000hz = model->L * 4000.0f / (SAMP_RATE / 2);

	float sig = 1e-4, elow = 1e-4, ehigh = 1e-4;
	for (int l = 1; l <= l_4000hz; l++)
	{
		float t = model->A[l] * model->A[l];

		if (l <= l_2000hz)
		{
			elow += t;
			if (l <= l_1000hz)
				sig += t;
		}
		if (l >= l_2000hz)
			ehigh += t;
	}

	float error = 1e-4; /* accumulated error between original and synthesised */

	/* Just test across the harmonics in the first 1000 Hz */
	for (int l = 1; l <= l_1000hz; l++)
	{
		complex_t Am; /* amplitude sample for this band */

		Am.r = 0.0f;
		Am.i = 0.0f;

		float den = 0.0f; /* denominator of Am expression */

		int al = ceilf_fast((l - 0.5f) * Wo_bin);
		int bl = ceilf_fast((l + 0.5f) * Wo_bin);

		/* Estimate amplitude of harmonic assuming harmonic is totally voiced */
		int offset = FFT_ENC / 2 - l * Wo_bin + 0.5; /* centers Hw[] about current harmonic */

		if (offset < -al)
			offset = -al;
		if (offset + bl >= FFT_ENC)
			offset = FFT_ENC - bl - 1;

		const float *restrict Wp = &W[offset];

		if (al < 0)
			al = 0;
		if (bl > FFT_ENC)
			bl = FFT_ENC;

		for (int m = al; m < bl; m++)
		{
			float w = Wp[m];

			Am.r += Sw[m].r * w;
			Am.i += Sw[m].i * w;
			den += w * w;
		}

		if (den < 1e-12f)
			continue;

		Am.r = Am.r / den;
		Am.i = Am.i / den;

		/* Determine error between estimated harmonic and original */
		for (int m = al; m < bl; m++)
		{
			float dr = Sw[m].r - Am.r * Wp[m];
			float di = Sw[m].i - Am.i * Wp[m];
			error += dr * dr + di * di;
		}
	}

	if (sig > V_THRESH_LIN * error)
		model->voiced = 1;
	else
		model->voiced = 0;

	/* post processing, helps clean up some voicing errors ------------------*/
	/* Look for Type 1 errors, strongly V speech that has been
	   accidentally declared UV */
	if (model->voiced == 0)
		if (elow > 10.0f * ehigh) // >10dB?
			model->voiced = 1;

	/* Look for Type 2 errors, strongly UV speech that has been
	   accidentally declared V */
	if (model->voiced == 1)
	{
		if (elow < 0.1f * ehigh) // <-10dB?
			model->voiced = 0;

		/* A common source of Type 2 errors is the pitch estimator
		   gives a low (50Hz) estimate for UV speech, which gives a
		   good match with noise due to the close harmoonic spacing.
		   These errors are much more common than people with 50Hz3
		   pitch, so we have just a small eratio threshold. */
		static const float sixty = 60.0f * TWO_PI / SAMP_RATE;
		if ((elow < 0.39810717055349725077f * ehigh) && (model->Wo <= sixty)) // <-4dB?
			model->voiced = 0;
	}
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
	make_analysis_window(c2, c2->fft_fwd_cfg, c2->w, c2->W);
}
