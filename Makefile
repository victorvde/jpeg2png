CC:=gcc
CFLAGS:=-std=c11 -pedantic -Ofast -msse2 -mfpmath=sse -DNDEBUG
WARN_FLAGS:=-Wall -Wextra -Winline
# CFLAGS:=-std=c11 -pedantic -Og -g
LFLAGS:= -s
BFLAGS:=
LIBS:=-ljpeg -lpng -lfftw3f -lm
OBJS:=jpeg2png.o utils.o jpeg.o png.o box.o upsample.o compute.o logger.o gopt/gopt.o

jpeg2png: $(OBJS)
	$(CC) $^ -o $@ $(LFLAGS) $(BFLAGS) $(LIBS)

-include $(OBJS:.o=.d)

gopt/gopt.o: gopt/gopt.c gopt/gopt.h
	$(CC) $< -c -o $@ $(CFLAGS) $(BFLAGS)

%.o: %.c
	$(CC) -MP -MMD $< -c -o $@ $(CFLAGS) $(BFLAGS) $(WARN_FLAGS)
