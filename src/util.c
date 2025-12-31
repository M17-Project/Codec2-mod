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

float fast_acosf(float x)
{
	/* Clamp to valid domain */
	if (x < -1.0f)
		x = -1.0f;
	if (x > 1.0f)
		x = 1.0f;

	int neg = 0;
	if (x < 0.0f)
	{
		x = -x;
		neg = 1;
	}

	/* Horner evaluation of polynomial */
	float p = ((((-0.004731162409361028f * x + 0.02013917916990988f) * x - 0.04547268716508862f) * x + 0.08798969877894365f) * x - 0.2145148619010314f) * x + 1.57079461300486f;
	float y = sqrtf(1.0f - x) * p;

	return neg ? (float)M_PI - y : y;
}

/* Note: the domain is restricted to [0, pi] */
/* this is enough for LSP space and avoids all the cosf() overhead */
float fast_cosf(float x)
{
	static const float PI_2 = M_PI / 2.0f;
	float sign = 1.0f;

	/* Fold to [0, pi/2] */
	if (x > PI_2)
	{				  /* pi/2 */
		x = M_PI - x; /* pi - x */
		sign = -1.0f;
	}

	/* Even polynomial Horner evaluation */
	float x2 = x * x;
	float p = (((2.323747198784698e-5f * x2 - 0.001385741491153336f) * x2 + 0.04166409053859771f) * x2 - 0.4999992686979083f) * x2 + 0.999999967262265f;

	return sign * p;
}

int codec2_rand(unsigned long *prng_state)
{
	*prng_state = *prng_state * 1103515245 + 12345;
	return ((unsigned)(*prng_state / 65536) % 32768);
}
