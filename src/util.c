#include "codec2_internal.h"
#include "util.h"

/* faster than atan2f(), but approximate (good enough) */
float fast_atan2f(float y, float x)
{
	static const float ONEQTR_PI = M_PI * 0.25f;
	static const float THRQTR_PI = M_PI * 0.75f;

	float r, angle;
	float abs_y = fabsf(y) + 1e-12f;

	if (x < 0.0f)
	{
		r = (x + abs_y) / (abs_y - x);
		angle = THRQTR_PI;
	}
	else
	{
		r = (x - abs_y) / (x + abs_y);
		angle = ONEQTR_PI;
	}
	angle += (0.1963f * r * r - 0.9817f) * r;
	return (y < 0.0f) ? -angle : angle;
}

int codec2_rand(unsigned long *prng_state)
{
	*prng_state = *prng_state * 1103515245 + 12345;
	return ((unsigned)(*prng_state / 65536) % 32768);
}
