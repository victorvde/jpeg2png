CC:=gcc
CFLAGS:=-std=c11 -Wall -Wextra -pedantic -O2 -ffast-math -march=pentium4 -ftree-vectorize -s
LIBS:=-ljpeg -lpng -lfftw3f -lm

jpeg2png: jpeg2png.c
	$(CC) $< -o $@ $(CFLAGS) $(LIBS)
