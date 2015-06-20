# Reset implicit rules as if using -r
.SUFFIXES:
# Reset implicit variables as if using -R
$(foreach var,$(filter-out .% MAKE% SUFFIXES,$(.VARIABLES)),$(if $(findstring $(origin $(var)),default),$(eval undefine $(var))))

# Build options
SSE2?=1
BUILTINS?=1
SIMD?=1
OPENMP?=1
DEBUG?=0
SAVE_ASM?=0
WINDOWS?=0

# VARIABLES
CFLAGS+=-std=c11 -pedantic
WARN_FLAGS+=-Wall -Wextra -Winline -Wshadow
CC?=$(HOST)gcc
WINDRES?=$(HOST)windres
LIBS+=-ljpeg -lpng -lm -lz
OBJS+=jpeg2png.o utils.o jpeg.o png.o box.o compute.o logger.o progressbar.o gopt/gopt.o ooura/dct.o

ifeq ($(SSE2),1)
CFLAGS+=-msse2 -mfpmath=sse
endif

ifeq ($(BUILTINS),1)
CFLAGS+=-DBUILTIN_UNREACHABLE -DBUILTIN_ASSUME_ALIGNED -DATTRIBUTE_UNUSED
endif

ifeq ($(SIMD),1)
CFLAGS+=-DUSE_SIMD
endif

ifeq ($(OPENMP),1)
CFLAGS+=-DUSE_OPENMP
BFLAGS+=-fopenmp
endif

ifeq ($(DEBUG),1)
CFLAGS+=-Og -DNDEBUG
BFLAGS+=-pg -g
else
CFLAGS+=-O3 -DNDEBUG
LDFLAGS+=-s
endif

ifeq ($(WINDOWS),1)
HOST?=i686-w64-mingw32-
EXE?=.exe
LDFLAGS+=-static

RES+=icon.rc.o
endif

ifeq ($(SAVE_ASM),1)
CFLAGS+=-save-temps -masm=intel -fverbose-asm
endif

# RULES
jpeg2png$(EXE): $(OBJS) $(RES) Makefile
	$(CC) $(OBJS) $(RES) -o $@ $(LDFLAGS) $(BFLAGS) $(LIBS)

-include $(OBJS:.o=.d)

gopt/gopt.o: gopt/gopt.c gopt/gopt.h Makefile
	$(CC) $< -c -o $@ $(CFLAGS) $(BFLAGS)

%.o: %.c Makefile
	$(CC) -MP -MMD $< -c -o $@ $(CFLAGS) $(BFLAGS) $(WARN_FLAGS)

%.rc.o: %.rc Makefile
	$(WINDRES) $< $@

.PHONY: clean
clean:
	git clean -Xf
