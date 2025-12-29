CC := gcc

CFLAGS_COMMON := -Iinc -Wall -Wextra -Wuninitialized -Wunused-parameter -std=c11 -D_GNU_SOURCE
CFLAGS_OPT := -O2
CFLAGS_MATH := #-ffast-math

CFLAGS := $(CFLAGS_COMMON) $(CFLAGS_OPT) $(CFLAGS_MATH)
LDFLAGS := -lm

MOD_SRC := src/*.c modified.c
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
