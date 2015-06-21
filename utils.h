#ifndef JPEG2PNG_UTILS_H
#define JPEG2PNG_UTILS_H

#include <assert.h>
#include <stdnoreturn.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

noreturn void die(const char *msg, ...);
noreturn void die_perror(const char *msg, ...);
clock_t start_timer(const char *name);
void stop_timer(clock_t t, const char *n);
void compare(const char *name, unsigned w, unsigned h, float *new, float *old);

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
#define DUMP(v, f) do { printf( #v " = " f "\n", v); } while(false)
#define DUMP_SIMD(r) do { __m128 _t = r; printf( #r " = %.9e,%.9e,%.9e,%.9e\n", _t[0], _t[1], _t[2], _t[3]); } while(false)
#define SWAP(type, x, y) do { type _t = x; x = y; y = _t; } while(false)

#if defined(NDEBUG) && defined(BUILTIN_UNREACHABLE)
  #define ASSUME(x) do { if(!(x)) { __builtin_unreachable(); } } while(false)
#else
  #define ASSUME(x) assert(x)
#endif
#if defined(NDEBUG) && defined(BUILTIN_ASSUME_ALIGNED)
  #define ASSUME_ALIGNED(x) x = __builtin_assume_aligned(x, 16)
#else
  #define ASSUME_ALIGNED(x) ASSUME((((uintptr_t)x) & 15) == 0)
#endif
#ifdef ATTRIBUTE_UNUSED
  #define POSSIBLY_UNUSED __attribute__((unused))
#else
  #define POSSIBLY_UNUSED
#endif

#define STRINGIFY(x) #x
#ifdef _OPENMP
#define OPENMP(x) _Pragma(STRINGIFY(omp x))
#else
#define OPENMP(x)
#endif

#define START_TIMER(n) clock_t macro_timer_##n = start_timer(#n);
#define STOP_TIMER(n) stop_timer(macro_timer_##n, #n);

inline void check(unsigned x, unsigned y, unsigned w, unsigned h) {
        ASSUME(x < w);
        ASSUME(y < h);
        (void) x;
        (void) y;
        (void) w;
        (void) h;
}

inline float *p(float *in, unsigned x, unsigned y, unsigned w, unsigned h) {
        check(x, y, w, h);
        return &in[y * w + x];
}

inline float sqr(float x) {
        return x * x;
}

inline float *alloc_real(size_t n) {
#if defined(_WIN32)
        float *f = _aligned_malloc(n * sizeof(float), 16);
#else
        float *f = aligned_alloc(16, n * sizeof(float));
#endif
        if(!f) { die("allocation error"); }
        return f;
}

inline void free_real(float *p) {
#ifdef _WIN32
        _aligned_free(p);
#else
        free(p);
#endif
}

#endif
