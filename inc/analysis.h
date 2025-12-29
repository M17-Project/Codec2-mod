#ifndef ANALYSIS_H
#define ANALYSIS_H

void analyse_one_frame(codec2_t *c2, model_t *model, const int16_t *speech);
void analysis_init(codec2_t *c2);

#endif /* ANALYSIS_H */