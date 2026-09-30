#ifndef CODEC2_MOD_UTIL_H
#define CODEC2_MOD_UTIL_H

#include <stdint.h>
#include <math.h>

float fast_atan2f(float y, float x);
float fast_acosf(float x);
float fast_cosf(float x);

#if defined(HAVE_SINCOSF)
#define codec2_sincosf(x, s, c) sincosf((x), (s), (c))
#else
static inline void codec2_sincosf(float x, float *s, float *c)
{
    *s = sinf(x);
    *c = cosf(x);
}
#endif

static inline int ceilf_fast(float x)
{
    int i = (int)x;
    return i + (i < x);
}

static inline int codec2_rand(uint32_t *prng_state)
{
	*prng_state = *prng_state * 1103515245U + 12345U;
	return ((unsigned)(*prng_state >> 16) & 0x7FFF);
}

/* x^y for x > 0 and |y*log2(x)| < 126
   ~2e-5 relative error, no libm calls, only mul/add + int operations */
static inline float fast_powf_pos(float x, float y)
{
	union { float f; uint32_t u; } v = {x};

	/* log2(x) = exponent + log2(1 + mantissa) */
	float e = (float)((int)(v.u >> 23) - 127);
	v.u = (v.u & 0x007FFFFFu) | 0x3F800000u;
	float m = v.f - 1.0f;
	float l2 = e + m * (1.441879776f + m * (-0.708863876f + m * (0.4152409602f + m * (-0.193510391f + m * 0.04526550458f))));

	/* 2^t = 2^floor(t) * 2^frac(t) */
	float t = y * l2;
	float fl = floorf(t);
	float f = t - fl;
	float p = 0.999999896f + f * (0.6931546159f + f * (0.2401407868f + f * (0.05586325839f + f * (0.008946224944f + f * 0.001895108641f))));
	v.f = p;
	v.u += (uint32_t)((int32_t)fl << 23);
	return v.f;
}

#endif /* CODEC2_MOD_UTIL_H */
