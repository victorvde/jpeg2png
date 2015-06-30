#ifdef DEBUG
  #ifdef _WIN32
    #error "enabling floating point exceptions for debugging is not supported on windows"
  #endif
  #define _GNU_SOURCE
  #include <fenv.h>
#endif

void enable_fp_exceptions() {
#ifdef DEBUG
        // glibc specific
        feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
#endif
}
