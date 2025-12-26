all:
	gcc -Wall -Wextra -Wuninitialized -O2 delta_lsp_cb.c kiss_fft.c kiss_fftr.c main.c -o test -lm
clean:
	rm test
