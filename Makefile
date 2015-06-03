CC:=gcc
CFLAGS:=-std=c11 -Wall -Wextra -Winline -pedantic -Ofast -march=pentium4 -DNDEBUG
#CFLAGS:=-std=c11 -Wall -Wextra -pedantic -Og -g
LFLAGS:=-s
LIBS:=-ljpeg -lpng -lfftw3f -lm
OBJS:=jpeg2png.o utils.o jpeg.o png.o box.o upsample.o compute.o

jpeg2png: $(OBJS)
	$(CC) $^ -o $@ $(LFLAGS) $(LIBS)

-include $(OBJS:.o=.d)

%.o: %.c
	$(CC) -MMD $< -c -o $@ $(CFLAGS)
