#ifndef JPEG2PNG_UTILS_H
#define JPEG2PNG_UTILS_H

#include <assert.h>
#include <stdnoreturn.h>
#include <time.h>
#include <math.h>

noreturn void die(const char *msg, ...);
noreturn void die_perror(const char *msg, ...);
clock_t start_timer(const char *name);
void stop_timer(clock_t t, const char *n);

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
#define DUMP(v, f) do { printf( #v " = " f "\n", v); } while(false)

#define START_TIMER(n) clock_t macro_timer_##n = start_timer(#n);
#define STOP_TIMER(n) stop_timer(macro_timer_##n, #n);

inline void check(int x, int y, int w, int h) {
        assert(0 <= x);
        assert(x < w);
        assert(0 <= y);
        assert(y < h);
        (void) x;
        (void) y;
        (void) w;
        (void) h;
}

inline float *p(float *in, int x, int y, int w, int h) {
        check(x, y, w, h);
        return &in[y * w + x];
}

inline float a(int n) {
        if(n == 0) {
                return 1./sqrt(2.);
        } else {
                return 1.;
        }
}

#endif
