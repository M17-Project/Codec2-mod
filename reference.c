#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <codec2.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(void)
{
	void *c2 = codec2_create(CODEC2_MODE_3200);

	uint8_t encoded[8] = {0};
	int16_t speech[160];

	for (uint8_t i=0; i< 160; i++)
		speech[i] = 0.5f * sinf(i/80.0f * 2.0f * M_PI);

	for (uint8_t j=0; j<10; j++)
	{
		codec2_encode(c2, encoded, speech);

		for (uint8_t i=0; i<8; i++)
			fprintf(stderr, "%02X ", encoded[i]);
		fprintf(stderr, "\n");
	}

	return 0;
}
