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

static uint8_t quantise(
	const float *cb, /* current VQ codebook */
	float *vec,		 /* vector to quantise */
	float *w,		 /* weighting vector */
	float *se		 /* accumulated squared error */
)
{
	float e; /* current error */
	float diff;

	float besti = 0;	/* best index so far */
	float beste = 1E32; /* best error so far */
	for (int j = 0; j < 32; j++)
	{
		e = 0.0;
		diff = cb[j] - vec[0];
		e = diff * w[0] * diff * w[0];
		if (e < beste)
		{
			beste = e;
			besti = j;
		}
	}

	*se += beste;

	return (besti);
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

void encode_lspds_scalar(uint16_t *indexes, float *lsp)
{
	float lsp_hz[LPC_ORD];
	float lsp__hz[LPC_ORD];
	float dlsp[LPC_ORD];
	float dlsp_[LPC_ORD];
	float wt[LPC_ORD];
	const float *cb;
	float se;

	for (int i = 0; i < LPC_ORD; i++)
	{
		wt[i] = 1.0;
	}

	/* convert from radians to Hz so we can use human readable
	   frequencies */
	for (int i = 0; i < LPC_ORD; i++)
		lsp_hz[i] = (4000.0 / M_PI) * lsp[i];

	wt[0] = 1.0;
	for (int i = 0; i < LPC_ORD; i++)
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

void decode_lspds_scalar(float *lsp_, const int *indexes)
{
	float lsp__hz[LPC_ORD];
	float dlsp_[LPC_ORD];
	const float *cb;

	for (int i = 0; i < LPC_ORD; i++)
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
