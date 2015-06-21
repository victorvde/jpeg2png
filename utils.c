#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"

noreturn void die(const char *msg, ...)  {
        fprintf(stderr, "jpeg2png: ");
        va_list l;
        va_start(l, msg);
        vfprintf(stderr, msg, l);
        va_end(l);
        fprintf(stderr, "\n");
        exit(EXIT_FAILURE);
}

noreturn void die_perror(const char *msg, ...)  {
        fprintf(stderr, "jpeg2png: ");
        va_list l;
        va_start(l, msg);
        vfprintf(stderr, msg, l);
        va_end(l);
        fprintf(stderr, ": ");
        perror(NULL);
        exit(EXIT_FAILURE);
}

clock_t start_timer(const char *name) {
        (void) name;
        return clock();
}

void stop_timer(clock_t t, const char *n) {
        clock_t diff = clock() - t;
        unsigned msec = (((double)diff / (double)CLOCKS_PER_SEC) * 1000.);
        printf("%s: %u ms\n", n, msec);
}

void compare(const char * name, unsigned w, unsigned h, float *new, float *old) {
        const float epsilon = 1.e-6;
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        float new1 = *p(new, x, y, w, h);
                        float old1 = *p(old, x, y, w, h);
                        if(isnan(new1)) {
                                die("%u, %u is NaN", x, y);
                        } else if (fabs(new1 - old1) / old1 > epsilon) {
                                die("difference at %s %u, %u: %.9e, %.9e\n", name, x, y, new1, old1);
                        }
                }
        }
}

static void check(unsigned x, unsigned y, unsigned w, unsigned h) {
        ASSUME(x < w);
        ASSUME(y < h);
        (void) x;
        (void) y;
        (void) w;
        (void) h;
}


float *p(float *in, unsigned x, unsigned y, unsigned w, unsigned h) {
        check(x, y, w, h);
        unsigned i = (y / 8) * (w / 8) + x / 8;
        unsigned j = 8 * (y % 8) + x % 8;
        return &in[i * 64 + j];
}
