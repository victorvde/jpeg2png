CC:=gcc
CFLAGS:=-std=c11 -Wall -Wextra -pedantic -O2 -s
LIBS:=-ljpeg -lpng

jpeg2png: jpeg2png.c
	$(CC) $< -o $@ $(LIBS)
