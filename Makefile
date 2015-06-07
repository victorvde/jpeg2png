CC?=gcc
EXE?=
CFLAGS+=-std=c11 -pedantic -msse2 -mfpmath=sse
WARN_FLAGS+=-Wall -Wextra -Winline
CFLAGS+=-DBUILTIN_UNREACHABLE -DBUILTIN_ASSUME_ALIGNED -DATTRIBUTE_UNUSED
CFLAGS+=-DUSE_SIMD # -DSIMD_VERIFY
CFLAGS+=-DUSE_OPENMP
CFLAGS+=-O3 -DNDEBUG
# CFLAGS+=-Og -DNDEBUG
# CFLAGS+=-save-temps -masm=intel -fverbose-asm
LFLAGS+=-s
# BFLAGS+=-pg -g
BFLAGS+=-fopenmp
LIBS:=-ljpeg -lpng -lm -lz
OBJS:=jpeg2png.o utils.o jpeg.o png.o box.o upsample.o compute.o logger.o progressbar.o gopt/gopt.o ooura/dct.o

jpeg2png$(EXE): $(OBJS)
	$(CC) $^ -o $@ $(LFLAGS) $(BFLAGS) $(LIBS)

-include $(OBJS:.o=.d)

gopt/gopt.o: gopt/gopt.c gopt/gopt.h
	$(CC) $< -c -o $@ $(CFLAGS) $(BFLAGS)

%.o: %.c
	$(CC) -MP -MMD $< -c -o $@ $(CFLAGS) $(BFLAGS) $(WARN_FLAGS)
