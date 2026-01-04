#ifndef CODEC2_MOD_QUANTISE_H
#define CODEC2_MOD_QUANTISE_H

void pack(uint8_t *bitArray, unsigned int *bitIndex,
          unsigned int field, unsigned int fieldWidth);

int unpack(const uint8_t *bitArray, unsigned int *bitIndex,
           unsigned int fieldWidth);

int encode_Wo(float Wo, uint8_t bits);
float decode_Wo(int index, int bits);
int encode_energy(float e, int bits);
float decode_energy(int index, int bits);
void encode_lspds_scalar(int *indexes, const float *lsp);
void decode_lspds_scalar(float *lsp_, const int *indexes);

#endif /* CODEC2_MOD_QUANTISE_H */
