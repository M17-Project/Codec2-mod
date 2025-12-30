#ifndef CODEC2_MOD_INTERNAL_H
#define CODEC2_MOD_INTERNAL_H

/*
  This work uses code written by David Rowe VK5DGR et al.
  https://github.com/drowe67/codec2
*/

#include <stdint.h>
#include <math.h>

#include "kiss_fft.h"
#include "kiss_fftr.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define TWO_PI (2.0f * M_PI)

#ifndef POW10F
#define POW10F(x) powf(10.0f, (x))
#endif

#define BYTES_PER_FRAME 8
#define SAMPLES_PER_FRAME (2 * N_SAMP)

#define N_S 0.01
#define MAX_AMP 80

#define FFT_ENC 512
#define FFT_DEC 512
#define V_THRESH 6.0
#define LPC_ORD 10

#define P_MAX_S 0.0200f

#define SAMP_RATE 8000
#define N_SAMP 80
#define M_PITCH 320
#define P_MIN 20
#define P_MAX 160
#define W0_MIN ((2.0 * M_PI) / 160.0)
#define W0_MAX ((2.0 * M_PI) / 20.0)
#define NW 279
#define TW 40

#define PMAX_M 320
#define COEFF 0.95
#define PE_FFT_SIZE 512
#define DEC 5
#define NDEC (M_PITCH / DEC)
#define T 0.1
#define CNLP 0.3
#define NLP_NTAP 48

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
#define LSP_DELTA1 0.01
#define BG_THRESH 40.0
#define BG_BETA 0.1
#define BG_MARGIN 6.0
#define LPC_ENERGY_FLOOR 1e-6f

#define FFT_R (TWO_PI / FFT_ENC)
#define FFT_1_R (FFT_ENC / TWO_PI)

#define CODEC2_RAND_MAX 32767

#define FFT_FWD_MEM_BYTES 4360
#define FFTR_MEM_BYTES 5408

/* asserts */
/* this is the exact Codec2 3200 bps configuration */
/* changing any of these parameters requires DSP math modification and buffer resizing */
_Static_assert(BYTES_PER_FRAME == 8, "Codec2 3200 only");
_Static_assert(LPC_ORD == 10, "Codec2 3200 assumes LPC_ORD == 10");
_Static_assert(MAX_AMP >= 80, "MAX_AMP must be >= 80 for Codec2 3200");
_Static_assert(FFT_ENC == 512, "Codec2 3200 assumes FFT_ENC == 512");
_Static_assert(FFT_DEC == FFT_ENC, "Encoder/decoder FFT sizes must match");
_Static_assert(PE_FFT_SIZE == 512, "Pitch estimator FFT size must be 512");
_Static_assert((M_PITCH % DEC) == 0, "M_PITCH must be divisible by DEC");
_Static_assert(NDEC == (M_PITCH / DEC), "NDEC must equal M_PITCH / DEC");
_Static_assert(NLP_NTAP == 48, "NLP FIR unrolling assumes NLP_NTAP == 48");
_Static_assert(WO_BITS <= 8, "encode_Wo() assumes <= 8-bit fields");
_Static_assert(E_BITS <= 8, "encode_energy() assumes <= 8-bit fields");

/* consts */
extern const float nlp_fir[NLP_NTAP];
extern const float bw_gamma[LPC_ORD + 1];
extern const float delta_lsp_cb[LPC_ORD][32];

typedef kiss_fft_cpx complex_t;

typedef struct model_t
{
    float Wo;
    int L;
    float A[MAX_AMP + 1];
    float phi[MAX_AMP + 1];
    int voiced;
} model_t;

typedef struct nlp_t
{
    float w[PMAX_M / DEC];
    float sq[PMAX_M];
    float mem_x, mem_y;
    float mem_fir[NLP_NTAP * 2];
    int mem_pos;
    kiss_fftr_cfg fftr_cfg;
    float fftr_buff[PE_FFT_SIZE];
    complex_t Fw[PE_FFT_SIZE / 2 + 1];
} nlp_t;

typedef struct codec2_t
{
    unsigned long next_rn;

    float w[M_PITCH];
    float W[FFT_ENC];
    float Pn[2 * N_SAMP];
    float Sn[M_PITCH];
    nlp_t nlp;

    float Sn_[2 * N_SAMP];
    float ex_phase;
    float bg_est;
    float prev_f0_enc;

    model_t prev_model_dec;
    float prev_lsps_dec[LPC_ORD];
    float prev_e_dec;

    kiss_fft_cfg fft_fwd_cfg;
    kiss_fftr_cfg fftr_fwd_cfg;
    kiss_fftr_cfg fftr_inv_cfg;
    kiss_fft_cfg phase_fft_fwd_cfg;
    kiss_fft_cfg phase_fft_inv_cfg;

    kiss_fft_cpx fft_buffer[FFT_ENC];

    /*
     * fft_buffer scratch usage:
     * - only one logical use at a time
     * - no overlapping lifetimes
     * - sizes are semantic (FFT_ENC, FFT_DEC)
     * Violating this will cause silent DSP corruption.
     */
    uint8_t fft_fwd_mem[FFT_FWD_MEM_BYTES];
    uint8_t fftr_fwd_mem[FFTR_MEM_BYTES];
    uint8_t fftr_inv_mem[FFTR_MEM_BYTES];
} codec2_t;

_Static_assert(sizeof(((codec2_t *)0)->fft_buffer) >= FFT_ENC * sizeof(kiss_fft_cpx), "fft_buffer too small for FFT_ENC scratch");

#endif /* CODEC2_MOD_INTERNAL_H */
