#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "utils.h"

noreturn void die(const char *msg, ...)  {
        if(msg) {
                fprintf(stderr, "jpeg2png: ");
                va_list l;
                va_start(l, msg);
                vfprintf(stderr, msg, l);
                va_end(l);
                fprintf(stderr, "\n");
        } else {
                perror("jpeg2png");
        }
        exit(EXIT_FAILURE);
}

clock_t start_timer(const char *name) {
        (void) name;
        return clock();
}

void stop_timer(clock_t t, const char *n) {
        clock_t diff = clock() - t;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("%s: %d ms\n", n, msec);
}
