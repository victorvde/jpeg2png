#ifndef JPEG2PNG_UTILS_H
#define JPEG2PNG_UTILS_H

#include <assert.h>
#include <stdnoreturn.h>
#include <time.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdint.h>

void die_message_start();
noreturn void die(const char *msg, ...);
noreturn void die_perror(const char *msg, ...);
clock_t start_timer(const char *name);
void stop_timer(clock_t t, const char *n);
void compare(const char *name, unsigned w, unsigned h, float *new, float *old);

// Convenience macros
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
#define DUMP(v, f) do { printf( #v " = " f "\n", v); } while(false)
#define DUMPD(v) DUMP(v, "%d")
#define DUMPF(v) DUMP(v, "%f")
#define DUMP_SIMD(r) do { __m128 _t = r; _mm_empty(); printf( #r " = %.9e,%.9e,%.9e,%.9e\n", _t[0], _t[1], _t[2], _t[3]); } while(false)
#define SWAP(type, x, y) do { type _t = x; x = y; y = _t; } while(false)
#define SQR(x) ((x) * (x))

// choose simd version or not, as chosen by the make parameter
#ifdef USE_SIMD
#define POSSIBLY_SIMD(x) x##_simd
#else
#define POSSIBLY_SIMD(x) x##_c
#endif

// asserts
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

// hide some warnings for some unused functions
#ifdef ATTRIBUTE_UNUSED
  #define POSSIBLY_UNUSED __attribute__((unused))
#else
  #define POSSIBLY_UNUSED
#endif

// nicer OpenMP pragmas
#define STRINGIFY(x) #x
#ifdef _OPENMP
#define OPENMP(x) _Pragma(STRINGIFY(omp x))
#else
#define OPENMP(x)
#endif

// timers for working on optimization
#define START_TIMER(n) clock_t macro_timer_##n = start_timer(#n);
#define STOP_TIMER(n) stop_timer(macro_timer_##n, #n);

// bounds check
static inline void check(unsigned x, unsigned y, unsigned w, unsigned h) {
        ASSUME(x < w);
        ASSUME(y < h);
        (void) x;
        (void) y;
        (void) w;
        (void) h;
}

// index image with bounds check
static inline float *p(float *in, unsigned x, unsigned y, unsigned w, unsigned h) {
        check(x, y, w, h);
        return &in[(size_t)y * w + x];
}

// convenience
static inline float sqf(float x) {
        return x * x;
}

// allocate aligned buffer for simd
static inline void *alloc_simd(size_t n) {
#if defined(_WIN32)
        void *p = _aligned_malloc(n, 16);
        if(!p) { die("allocation error"); }
#else
        void *p = NULL;
        if (posix_memalign(&p, 16, n) != 0) {
            die("aligned allocation error");
        }
#endif
        ASSUME_ALIGNED(p);
        return p;
}

static inline void free_simd(void *p) {
#ifdef _WIN32
        _aligned_free(p);
#else
        free(p);
#endif
}

#endif
