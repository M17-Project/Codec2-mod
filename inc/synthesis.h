#ifndef SYNTHESIS_H
#define SYNTHESIS_H

void synthesise_one_frame(
    codec2_t *c2,
    int16_t *speech,
    model_t *model,
    const complex_t *Aw,
    float gain
);

void synthesis_init(codec2_t *c2);

#endif /* SYNTHESIS_H */