#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <codec2.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
	struct timespec tick, tock;

	struct CODEC2 *c2;
	c2 = codec2_create(CODEC2_MODE_3200);

	uint8_t encoded[8] = {0};
	int16_t speech[160];

	FILE *audio_in, *bitstream, *audio_out;

	audio_in = fopen("../../sample.raw", "rb");
	bitstream = fopen("../../reference_bitstream.bin", "wb");
	audio_out = fopen("../../reference_decoded.raw", "wb");

	clock_gettime(CLOCK_MONOTONIC, &tick);
	size_t rb;
	while ((rb = fread(speech, 1, sizeof(speech), audio_in)) > 0)
	{
		if (rb < sizeof(speech))
			memset(&speech[rb / 2], 0, sizeof(speech) - rb); //`break;` here instead to cut away the remaining audio (just like `c2enc` does)

		codec2_encode(c2, encoded, speech);
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

	bitstream = fopen("../../reference_bitstream.bin", "rb");
	codec2_destroy(c2);
	c2 = codec2_create(CODEC2_MODE_3200);

	while (fread(encoded, 8, 1, bitstream) == 1)
	{
		codec2_decode(c2, speech, encoded);
		fwrite(speech, sizeof(speech), 1, audio_out);
	}

	fclose(audio_in);
	fclose(bitstream);
	fclose(audio_out);

	/*for (uint8_t i=0; i< 160; i++)
		speech[i] = 0.5f * sinf(i/80.0f * 2.0f * M_PI);

	for (uint8_t j=0; j<10; j++)
	{
		codec2_encode(c2, encoded, speech);

		for (uint8_t i=0; i<8; i++)
			fprintf(stderr, "%02X ", encoded[i]);
		fprintf(stderr, "\n");
	}*/

	return 0;
}
