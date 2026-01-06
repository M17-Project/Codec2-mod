#include "codec2_internal.h"
#include "quantise.h"

static inline unsigned int binary_to_gray(unsigned int x)
{
	return x ^ (x >> 1);
}

static inline unsigned int gray_to_binary(unsigned int g)
{
	g ^= g >> 16;
	g ^= g >> 8;
	g ^= g >> 4;
	g ^= g >> 2;
	g ^= g >> 1;
	return g;
}

void pack(
	uint8_t *bitArray,		/* The output bit array. */
	unsigned int *bitIndex, /* Index into the string in BITS, not bytes. */
	unsigned int field,		/* The bit field to be packed. */
	unsigned int fieldWidth /* Width of the field in BITS, not bytes. */
)
{
	/* Convert to Gray code if multi-bit */
	if (fieldWidth > 1)
		field = binary_to_gray(field);

	do
	{
		unsigned int bI = *bitIndex;
		unsigned int bitsLeft = 8 - (bI & 0x7);
		unsigned int sliceWidth = bitsLeft < fieldWidth ? bitsLeft : fieldWidth;
		unsigned int wordIndex = bI >> 3;

		bitArray[wordIndex] |=
			(uint8_t)((field >> (fieldWidth - sliceWidth)) << (bitsLeft - sliceWidth));

		*bitIndex = bI + sliceWidth;
		fieldWidth -= sliceWidth;

	} while (fieldWidth != 0);
}

int unpack(
	const uint8_t *bitArray, /* The input bit string. */
	unsigned int *bitIndex,	 /* Index into the string in BITS, not bytes. */
	unsigned int fieldWidth	 /* Width of the field in BITS, not bytes. */
)
{
	static const uint8_t mask8[9] = {
		0x00, 0x01, 0x03, 0x07, 0x0F,
		0x1F, 0x3F, 0x7F, 0xFF};

	unsigned int field = 0;
	unsigned int origWidth = fieldWidth;

	do
	{
		unsigned int bI = *bitIndex;
		unsigned int bitsLeft = 8 - (bI & 0x7);
		unsigned int sliceWidth = bitsLeft < fieldWidth ? bitsLeft : fieldWidth;

		field |= (((bitArray[bI >> 3] >> (bitsLeft - sliceWidth)) &
				   mask8[sliceWidth])
				  << (fieldWidth - sliceWidth));

		*bitIndex = bI + sliceWidth;
		fieldWidth -= sliceWidth;
	} while (fieldWidth != 0);

	return (origWidth > 1) ? gray_to_binary(field) : field;
}

int encode_Wo(float Wo, uint8_t bits)
{
	int Wo_levels = 1 << bits;

	float norm = (Wo - W0_MIN) / (W0_MAX - W0_MIN);
	int index = (int)(Wo_levels * norm + 0.5f);

	if (index < 0)
		index = 0;
	else if (index > Wo_levels - 1)
		index = Wo_levels - 1;

	return index;
}

float decode_Wo(int index, int bits)
{
	int Wo_levels = 1 << bits;
	float step = (W0_MAX - W0_MIN) / Wo_levels;

	// note: the return value should be
	// `W0_MIN + step * (index + 0.5f)`
	// to achieve symmetry with the encoder,
	// but for some reason, it lowers ViSQOL MOS slightly
	return W0_MIN + step * index;
}

int encode_energy(float e, int bits)
{
	const int e_levels = 1 << bits;
	const float e_min = E_MIN_DB;
	const float e_max = E_MAX_DB;

	float e_db = 10.0f * log10f(fmaxf(e, 1e-12f)); // avoid log(0)
	float norm = (e_db - e_min) / (e_max - e_min);

	int index = (int)(e_levels * norm + 0.5f);

	if (index < 0)
		index = 0;
	else if (index > e_levels - 1)
		index = e_levels - 1;

	return index;
}

float decode_energy(int index, int bits)
{
	const int e_levels = 1 << bits;
	const float step = (E_MAX_DB - E_MIN_DB) / e_levels;

	// note: just like decode_Wo() above
	// (index + 0.5f) would assure symmetry with the encoder
	// but for some reason it breaks MOS
	float e_db = E_MIN_DB + step * index;
	return POW10F(e_db / 10.0f);
}

void encode_lspds_scalar(int *indexes, const float *lsp)
{
	static const float k = 4000.0f / M_PI;
	int last_q_hz = 0;

	for (int i = 0; i < LPC_ORD; i++)
	{
		const uint16_t *cb = delta_lsp_cb[i];
		float lsp_hz = k * lsp[i];
		float target = (i == 0) ? lsp_hz : (lsp_hz - (float)last_q_hz);

		// given the Poisson-like probability density function of the LSF indices
		// and strictly monotonic codebooks,
		// it seems reasonable to use binary search
		// indices for LSF 7, 8, and 9 (counting from 1)
		// have high peaks at cb[31], so we can detect that early
		// modes for LSF indices: {8, 4, 6, 11, 10, 8, 15, 13, 14, 11}
		// (measured using 500 random files from LibriSpeech)
		int best_j;

		/* early saturation detection */
		if (target <= (float)cb[0])
		{
			best_j = 0;
		}
		else if (target >= (float)cb[31])
		{
			best_j = 31;
		}
		else
		{
			/* find first j such that cb[j] > target */
			int lo = 0, hi = 31;
			while (lo + 1 < hi)
			{
				int mid = (lo + hi) >> 1;
				if ((float)cb[mid] <= target)
					lo = mid;
				else
					hi = mid;
			}

			float e_lo = fabsf((float)cb[lo] - target);
			float e_hi = fabsf((float)cb[hi] - target);

			/* lower index wins ties */
			best_j = (e_lo <= e_hi) ? lo : hi;
		}

		indexes[i] = best_j;

		if (i == 0)
			last_q_hz = cb[best_j];
		else
			last_q_hz += cb[best_j];
	}
}

void decode_lspds_scalar(float *lsp_, const int *indexes)
{
	int lsp_hz = 0;
	static const float k = M_PI / 4000.0f;

	for (int i = 0; i < LPC_ORD; i++)
	{
		lsp_hz += delta_lsp_cb[i][indexes[i]];
		lsp_[i] = k * (float)lsp_hz;
	}
}
