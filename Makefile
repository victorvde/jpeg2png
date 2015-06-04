CC:=gcc
CFLAGS:=-std=c11 -pedantic -Ofast -march=pentium4 -DNDEBUG
WARN_FLAGS:=-Wall -Wextra -Winline
#CFLAGS:=-std=c11 -Wall -Wextra -pedantic -Og -g
LFLAGS:=-s
LIBS:=-ljpeg -lpng -lfftw3f -lm
OBJS:=jpeg2png.o utils.o jpeg.o png.o box.o upsample.o compute.o logger.o gopt/gopt.o

jpeg2png: $(OBJS)
	$(CC) $^ -o $@ $(LFLAGS) $(LIBS)

-include $(OBJS:.o=.d)

gopt/gopt.o: gopt/gopt.c gopt/gopt.h
	$(CC) $< -c -o $@ $(CFLAGS)

%.o: %.c
	$(CC) -MP -MMD $< -c -o $@ $(CFLAGS) $(WARN_FLAGS)
