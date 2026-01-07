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

#endif /* CODEC2_MOD_UTIL_H */
