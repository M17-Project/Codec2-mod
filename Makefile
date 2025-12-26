all:
	gcc -Wall -Wextra -Wuninitialized -O2 main.c kiss_fft.c kiss_fftr.c -o test -lm
clean:
	rm *.o test
