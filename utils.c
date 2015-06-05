#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

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
        unsigned msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("%s: %u ms\n", n, msec);
}
