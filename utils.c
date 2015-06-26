#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "progressbar.h"
#include "jpeg2png.h"

static void die_common(const char *msg, va_list l) {
        if(main_progressbar) {
                progressbar_done(main_progressbar);
                main_progressbar = NULL;
        }
        fprintf(stderr, "jpeg2png: ");
        vfprintf(stderr, msg, l);
}

noreturn void die(const char *msg, ...)  {
        va_list l;
        va_start(l, msg);
        die_common(msg, l);
        va_end(l);
        fprintf(stderr, "\n");
        exit(EXIT_FAILURE);
}

noreturn void die_perror(const char *msg, ...)  {
        va_list l;
        va_start(l, msg);
        die_common(msg, l);
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
