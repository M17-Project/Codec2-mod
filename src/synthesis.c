#include "codec2_internal.h"
#include "synthesis.h"
#include "util.h"
#include <string.h>
#include <math.h>

static void make_synthesis_window(float *restrict Pn)
{
	static const int n0 = N_SAMP / 2;
	static const int n1 = 3 * N_SAMP / 2;
	static const float inv_2TW = 1.0f / (2.0f * TW);

	/* Generate Parzen window in time domain */
	/* zero */
	memset(Pn, 0, (n0 - TW) * sizeof(float));

	/* ramp up */
	for (int i = n0 - TW; i < n0 + TW; i++)
		Pn[i] = (i - (n0 - TW)) * inv_2TW;

	/* flat */
	for (int i = n0 + TW; i < n1 - TW; i++)
		Pn[i] = 1.0f;

	/* ramp down */
	for (int i = n1 - TW; i < n1 + TW; i++)
		Pn[i] = ((n1 + TW) - i) * inv_2TW;

	/* zero */
	memset(&Pn[n1 + TW], 0, (2 * N_SAMP - (n1 + TW)) * sizeof(float));
}

static void sample_phase(
	const model_t *restrict model,
	complex_t *restrict H,
	const complex_t *restrict A /* LPC analysis filter in freq domain */
)
{
	const float k = model->Wo / FFT_R; // pre-compute

	/* Sample phase at harmonics */
	for (int m = 1; m <= model->L; m++)
	{
		int b = (int)(m * k + 0.5);

		// clamp b
		if (b >= FFT_ENC / 2)
			b = FFT_ENC / 2 - 1;

		/* synth filter 1/A is opposite phase to analysis filter */
		// H[m] = cconj(A[b]);
		H[m].r = A[b].r;
		H[m].i = -A[b].i;
	}
}

static void phase_synth_zero_order(
	codec2_t *c2,
	model_t *model,
	float *ex_phase,   /* excitation phase of fundamental        */
	const complex_t *H /* L synthesis filter freq domain samples */
)
{
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
	float phi0 = ex_phase[0];

	static const float k = TWO_PI / CODEC2_RAND_MAX;

	for (int m = 1; m <= model->L; m++)
	{
		/* generate excitation */
		if (model->voiced)
		{
			float s, c;
			codec2_sincosf(phi0 * m, &s, &c);
			Ex[m].r = c;
			Ex[m].i = s;
		}
		else
		{
			/* When a few samples were tested I found that LPC filter
			   phase is not needed in the unvoiced case, but no harm in
			   keeping it.
			*/
			float phi = k * (float)codec2_rand(&c2->next_rn);
			float s, c;
			codec2_sincosf(phi, &s, &c);
			Ex[m].r = c;
			Ex[m].i = s;
		}

		/* filter using LPC filter */
		A_[m].r = H[m].r * Ex[m].r - H[m].i * Ex[m].i;
		A_[m].i = H[m].i * Ex[m].r + H[m].r * Ex[m].i;

		/* modify sinusoidal phase */
		new_phi = fast_atan2f(A_[m].i, A_[m].r + 1e-12);
		model->phi[m] = new_phi;
	}
}

static void postfilter(codec2_t *restrict c2, model_t *restrict model, float *restrict bg_est)
{
	static const float k = TWO_PI / CODEC2_RAND_MAX;

	/* determine average energy across spectrum */
	float e = 1e-12;
	for (int m = 1; m <= model->L; m++)
		e += model->A[m] * model->A[m];

	e = 10.0 * log10f(e / model->L);

	/* If beneath threshold, update bg estimate.  The idea
	   of the threshold is to prevent updating during high level
	   speech. */
	if (!model->voiced)
	{
		if (e < BG_THRESH)
			*bg_est = *bg_est * (1.0 - BG_BETA) + e * BG_BETA;
	}
	else
	{
		/* now mess with phases during voiced frames to make any harmonics
		   less then our background estimate unvoiced. */

		float thresh = POW10F((*bg_est + BG_MARGIN) / 20.0f);

		for (int m = 1; m <= model->L; m++)
		{
			if (model->A[m] < thresh)
			{
				model->phi[m] = k * (float)codec2_rand(&c2->next_rn);
			}
		}
	}
}

static void ear_protection(float *in_out, int n)
{
	float max_abs = 0.0f;

	/* find maximum sample in frame */
	for (int i = 0; i < n; i++)
	{
		float tmp = fabsf(in_out[i]);
		if (tmp > max_abs)
			max_abs = tmp;
	}

	if (max_abs <= 30000.0f)
		return; // nothing to do here

	/* determine how far above set point */
	float over = max_abs / 30000.0f;

	/* If we are x dB over set point we reduce level by 2x dB, this
	   attenuates major excursions in amplitude (likely to be caused
	   by bit errors) more than smaller ones */
	if (over > 1.0f)
	{
		float gain = 1.0f / (over * over);
		for (int i = 0; i < n; i++)
			in_out[i] *= gain;
	}
}

static void synthesise(
	codec2_t *c2,
	kiss_fftr_cfg fftr_inv_cfg,
	float *Sn_,					   /* time domain synthesised signal              */
	const model_t *restrict model, /* ptr to model parameters for this frame      */
	const float *restrict Pn,	   /* time domain Parzen window                   */
	int shift					   /* flag used to handle transition frames       */
)
{
	// NOTE: lifetimes do not overlap
	complex_t *Sw_ = c2->fft_buffer;	  /* DFT of synthesised signal */
	float *sw_ = (float *)c2->fft_buffer; /* synthesised signal */

	if (shift)
	{
		/* Update memories */
		memmove(Sn_, &Sn_[N_SAMP], (N_SAMP - 1) * sizeof(float));
		Sn_[N_SAMP - 1] = 0.0;
	}

	memset(Sw_, 0, (FFT_DEC / 2 + 1) * sizeof(complex_t)); // original Sw_ size was this

	/* Now set up frequency domain synthesised speech */
	const float Wo_bin = model->Wo * FFT_1_R;

	for (int l = 1; l <= model->L; l++)
	{
		int b = (int)(l * Wo_bin + 0.5); // FFT_DEC == FFT_ENC
		if (b > ((FFT_DEC / 2) - 1))
		{
			b = (FFT_DEC / 2) - 1;
		}
		float s, c;
		codec2_sincosf(model->phi[l], &s, &c);
		Sw_[b].r = model->A[l] * c;
		Sw_[b].i = model->A[l] * s;
	}

	/* Perform inverse DFT */
	kiss_fftri(fftr_inv_cfg, Sw_, sw_);

	/* Overlap add to previous samples */
	for (int i = 0; i < N_SAMP - 1; i++)
	{
		Sn_[i] += sw_[FFT_DEC - N_SAMP + 1 + i] * Pn[i];
	}

	float *restrict dst = &Sn_[N_SAMP - 1];
	const float *restrict src = sw_;
	const float *restrict win = &Pn[N_SAMP - 1];

	if (shift)
	{
		for (int j = 0; j < N_SAMP + 1; j++)
			dst[j] = src[j] * win[j];
	}
	else
	{
		for (int j = 0; j < N_SAMP + 1; j++)
			dst[j] += src[j] * win[j];
	}
}

void synthesise_one_frame(
	codec2_t *c2,
	int16_t *speech,
	model_t *model,
	const complex_t *Aw,
	float gain)
{
	/* LPC based phase synthesis */
	complex_t *H = (complex_t *)c2->fft_buffer; // use a chunk of that array as scratch
	sample_phase(model, H, Aw);
	phase_synth_zero_order(c2, model, &c2->ex_phase, H);
	postfilter(c2, model, &c2->bg_est);
	synthesise(c2, c2->fftr_inv_cfg, c2->Sn_, model, c2->Pn, 1);

	for (int i = 0; i < N_SAMP; i++)
	{
		c2->Sn_[i] *= gain;
	}

	ear_protection(c2->Sn_, N_SAMP);

	for (int i = 0; i < N_SAMP; i++)
	{
		if (c2->Sn_[i] > 32767.0f)
			speech[i] = 32767;
		else if (c2->Sn_[i] < -32767.0f)
			speech[i] = -32767;
		else
			speech[i] = c2->Sn_[i];
	}
}

void synthesis_init(codec2_t *c2)
{
	make_synthesis_window(c2->Pn);
}
