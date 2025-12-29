#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <time.h>

#include "codec2_mod.h"

int main(void)
{
	struct timespec tick, tock;

	codec2_t c2;
	codec2_init(&c2);

	uint8_t encoded[CODEC2_BYTES_PER_FRAME];
	int16_t speech[CODEC2_SAMPLES_PER_FRAME];

	FILE *audio_in, *bitstream, *audio_out;

	audio_in = fopen("../../sample.raw", "rb");
	bitstream = fopen("../../modified_bitstream.bin", "wb");
	audio_out = fopen("../../modified_decoded.raw", "wb");

	clock_gettime(CLOCK_MONOTONIC, &tick);
	size_t rb;
	while ((rb = fread(speech, 1, sizeof(speech), audio_in)) > 0)
	{
		if (rb < sizeof(speech))
			memset(&speech[rb / 2], 0, sizeof(speech) - rb); //`break;` here instead to cut away the remaining audio (just like `c2enc` does)

		codec2_encode(&c2, encoded, speech);
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

	codec2_init(&c2);
	bitstream = fopen("../../modified_bitstream.bin", "rb");

	while (fread(encoded, CODEC2_BYTES_PER_FRAME, 1, bitstream) == 1)
	{
		codec2_decode(&c2, speech, encoded);
		fwrite(speech, sizeof(speech), 1, audio_out);
	}

	fclose(audio_in);
	fclose(bitstream);
	fclose(audio_out);

	/*for (uint8_t i = 0; i < 160; i++)
		speech[i] = 0.5f * sinf(i / 80.0f * TWO_PI);

	for (uint8_t j = 0; j < 10; j++)
	{
		codec2_encode(&c2, encoded, speech);

		for (uint8_t i = 0; i < 8; i++)
			fprintf(stderr, "%02X ", encoded[i]);
		fprintf(stderr, "\n");
	}*/

	return 0;
}
