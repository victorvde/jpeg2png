CC:=gcc
CFLAGS:=-std=c11 -Wall -Wextra -Winline -pedantic -Ofast -march=pentium4 -s -DNDEBUG
#CFLAGS:=-std=c11 -Wall -Wextra -pedantic -Og -g
LIBS:=-ljpeg -lpng -lfftw3f -lm

SRCS:=jpeg2png.c utils.c jpeg.c

jpeg2png: jpeg2png.o utils.o jpeg.o png.o
	$(CC) $^ -o $@ $(LIBS)

%o: %.c *.h
	$(CC) $< -o $@ $(CFLAGS)
