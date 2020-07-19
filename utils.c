#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "progressbar.h"
#include "jpeg2png.h"

// clean up line and print prefix
void die_message_start() {
        if(main_progressbar) {
                progressbar_clear(main_progressbar);
                main_progressbar = NULL;
        }
        fprintf(stderr, "jpeg2png: ");
}

// abort program with a message
noreturn void die(const char *msg, ...)  {
        die_message_start();
        va_list l;
        va_start(l, msg);
        vfprintf(stderr, msg, l);
        va_end(l);
        fprintf(stderr, "\n");
        exit(EXIT_FAILURE);
}

// abort program with a system message
noreturn void die_perror(const char *msg, ...)  {
        die_message_start();
        va_list l;
        va_start(l, msg);
        vfprintf(stderr, msg, l);
        va_end(l);
        fprintf(stderr, ": ");
        perror(NULL);
        exit(EXIT_FAILURE);
}

// see utils.h
clock_t start_timer(const char *name) {
        (void) name;
        return clock();
}

void stop_timer(clock_t t, const char *n) {
        clock_t diff = clock() - t;
        unsigned msec = (((double)diff / (double)CLOCKS_PER_SEC) * 1000.);
        printf("%s: %u ms\n", n, msec);
}

// compare image sized buffers, e.g. c and simd versions
void compare(const char * name, unsigned w, unsigned h, float *bnew, float *bold) {
        const float epsilon = 1.e-6;
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        float bnew1 = *p(bnew, x, y, w, h);
                        float bold1 = *p(bold, x, y, w, h);
                        if(isnan(bnew1)) {
                                die("%u, %u is NaN", x, y);
                        } else if (fabsf(bnew1 - bold1) / bold1 > epsilon) {
                                die("difference at %s %u, %u: %.9e, %.9e\n", name, x, y, bnew1, bold1);
                        }
                }
        }
}
