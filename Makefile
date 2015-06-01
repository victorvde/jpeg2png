CC:=gcc
CFLAGS:=-std=c11 -Wall -Wextra -Winline -pedantic -Ofast -march=pentium4 -s -DNDEBUG
#CFLAGS:=-std=c11 -Wall -Wextra -pedantic -Og -g
LIBS:=-ljpeg -lpng -lfftw3f -lm

jpeg2png: jpeg2png.c
	$(CC) $< -o $@ $(CFLAGS) $(LIBS)
