CC:=gcc
CFLAGS:=-std=c11 -Wall -Wextra -pedantic -O3 -ffast-math -s
LIBS:=-ljpeg -lpng -lm

jpeg2png: jpeg2png.c
	$(CC) $< -o $@ $(LIBS)
