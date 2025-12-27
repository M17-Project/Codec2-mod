CC := gcc

CFLAGS_COMMON := -Wall -Wextra -Wuninitialized -Wunused-parameter -std=c11 -D_GNU_SOURCE
CFLAGS_FP := -fno-fast-math -ffp-contract=off
CFLAGS_OPT := -O2
CFLAGS := $(CFLAGS_COMMON) $(CFLAGS_FP) $(CFLAGS_OPT)

LDFLAGS := -lm

MOD_SRC := delta_lsp_cb.c kiss_fft.c kiss_fftr.c modified.c
REF_SRC := reference.c

MOD_BIN := modified
REF_BIN := reference

all: $(MOD_BIN) $(REF_BIN)

$(MOD_BIN): $(MOD_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS)

$(REF_BIN): $(REF_SRC)
	$(CC) $(CFLAGS) $^ -o $@ $(LDFLAGS) -lcodec2

clean:
	rm -f $(MOD_BIN) $(REF_BIN)

.PHONY: all clean
