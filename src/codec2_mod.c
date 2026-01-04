#include "codec2_mod.h"
#include "codec2_internal.h"
#include "analysis.h"
#include "synthesis.h"
#include "nlp.h"
#include "lpc.h"
#include "interp.h"
#include "quantise.h"

void codec2_init(codec2_t *c2)
{
	c2->next_rn = 1; // random number geterator - seed

	for (int i = 0; i < M_PITCH; i++)
		c2->Sn[i] = 1.0f;

	memset(c2->Sn_, 0, sizeof(c2->Sn_));

	/* static FFT mem allocations */
	size_t mem;

	/* FFT forward */
	mem = sizeof(c2->fft_fwd_mem);
	c2->fft_fwd_cfg = kiss_fft_alloc(FFT_ENC, 0, c2->fft_fwd_mem, &mem);

	/* FFT real forward */
	mem = sizeof(c2->fftr_fwd_mem);
	c2->fftr_fwd_cfg = kiss_fftr_alloc(FFT_ENC, 0, c2->fftr_fwd_mem, &mem);

	/* FFT real inverse */
	mem = sizeof(c2->fftr_inv_mem);
	c2->fftr_inv_cfg = kiss_fftr_alloc(FFT_DEC, 1, c2->fftr_inv_mem, &mem);

	/* NLP FFT - reused (same type, direction, and size) */
	c2->nlp.fftr_cfg = c2->fftr_fwd_cfg;

	analysis_init(c2);
	synthesis_init(c2);

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
}

void codec2_encode(codec2_t *c2, uint8_t *bits, const int16_t *speech)
{
	model_t model = {0};
	float ak[LPC_ORD + 1];
	float lsps[LPC_ORD];
	float e;
	int Wo_index, e_index;
	int lspd_indexes[LPC_ORD];
	unsigned int nbit = 0;

	memset(bits, 0, BYTES_PER_FRAME); // 64 bits

	/* first 10ms analysis frame - we just want voicing */
	analyse_one_frame(c2, &model, speech);
	pack(bits, &nbit, model.voiced, 1);

	/* second 10ms analysis frame */
	analyse_one_frame(c2, &model, &speech[N_SAMP]);
	pack(bits, &nbit, model.voiced, 1);
	Wo_index = encode_Wo(model.Wo, WO_BITS);
	pack(bits, &nbit, Wo_index, WO_BITS);

	speech_to_uq_lsps(c2, lsps, ak, &e, c2->Sn, c2->w);
	e_index = encode_energy(e, E_BITS);
	pack(bits, &nbit, e_index, E_BITS);

	encode_lspds_scalar(lspd_indexes, lsps);
	for (int i = 0; i < LSPD_SCALAR_INDEXES; i++)
	{
		pack(bits, &nbit, lspd_indexes[i], 5); // codebook indices are 5 bits
	}
}

void codec2_decode(codec2_t *c2, int16_t *speech, const uint8_t *bits)
{
	model_t model[2] = {0};
	int lspd_indexes[LPC_ORD];
	float lsps[2][LPC_ORD];
	int Wo_index, e_index;
	float e[2];
	float ak[2][LPC_ORD + 1];
	complex_t Aw[FFT_ENC];
	float A2[FFT_ENC / 2];

	/* unpack bits from channel ------------------------------------*/
	unsigned int nbit = 0;
	/* this will partially fill the model params for the 2 x 10ms
	   frames */
	model[0].voiced = unpack(bits, &nbit, 1);
	model[1].voiced = unpack(bits, &nbit, 1);

	Wo_index = unpack(bits, &nbit, WO_BITS);
	model[1].Wo = decode_Wo(Wo_index, WO_BITS);
	model[1].L = M_PI / model[1].Wo;

	e_index = unpack(bits, &nbit, E_BITS);
	e[1] = decode_energy(e_index, E_BITS);

	for (int i = 0; i < LSPD_SCALAR_INDEXES; i++)
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

	for (int i = 0; i < 2; i++)
	{
		lsp_to_lpc(&lsps[i][0], &ak[i][0]);
		aks_to_mag2(c2, &ak[i][0], &model[i], e[i], Aw, A2);
		apply_lpc_correction(&model[i]);
		synthesise_one_frame(c2, &speech[N_SAMP * i], &model[i], Aw, 1.0f);
	}

	/* update memories for next frame ----------------------------*/
	c2->prev_model_dec = model[1];
	c2->prev_e_dec = e[1];
	for (int i = 0; i < LPC_ORD; i++)
		c2->prev_lsps_dec[i] = lsps[1][i];
}
