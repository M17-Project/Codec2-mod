#include "codec2_internal.h"
#include "quantise.h"

void pack(
	uint8_t *bitArray,		/* The output bit array. */
	unsigned int *bitIndex, /* Index into the string in BITS, not bytes. */
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
	unsigned int *bitIndex,	 /* Index into the string in BITS, not bytes. */
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

void encode_lspds_scalar(int *indexes, const float *lsp)
{
	float last_q_hz = 0.0f;

	for (int i = 0; i < LPC_ORD; i++)
	{
		const float *cb = delta_lsp_cb[i];
		float lsp_hz = (4000.0f / M_PI) * lsp[i];
		float target = (i == 0) ? lsp_hz : (lsp_hz - last_q_hz);

		float best_e = 4000.0f; // the value can't be more than Nyquist
		int best_j = 0;

		for (int j = 0; j < 32; j++)
		{
			float diff = fabsf(cb[j] - target);

			if (diff < best_e)
			{
				best_e = diff;
				best_j = j;
			}
			else if (cb[j] > target)
			{
				break;
			}
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
	float lsp__hz[LPC_ORD];
	float dlsp_[LPC_ORD];
	const float *cb;

	for (int i = 0; i < LPC_ORD; i++)
	{
		cb = delta_lsp_cb[i];
		dlsp_[i] = cb[indexes[i]];

		if (i)
			lsp__hz[i] = lsp__hz[i - 1] + dlsp_[i];
		else
			lsp__hz[0] = dlsp_[0];

		lsp_[i] = (M_PI / 4000.0f) * lsp__hz[i];
	}
}
