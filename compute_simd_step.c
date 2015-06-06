#include <xmmintrin.h>
#include <math.h>
#include "utils.h"

static double compute_step_tv_simd(unsigned w, unsigned h, float *in, float *objective_gradient, float *in_x, float *in_y)  {
        ASSUME_ALIGNED(in);
        ASSUME_ALIGNED(objective_gradient);
        ASSUME_ALIGNED(in_x);
        ASSUME_ALIGNED(in_y);

        ASSUME(w % 8 == 0);
        ASSUME(h % 8 == 0);

        unsigned simd_w = w - 4;

        double tv = 0.;

        const __m128 mm_inf = _mm_set_ps1(INFINITY);
        const __m128 mm_zero = _mm_set_ps1(0.);

        for(unsigned y = 0; y < h-1; y++) {
                unsigned y_w = y * w;
                for(unsigned x = 0; x < simd_w; x+=4) {
                        // load current pixels
                        __m128 here = _mm_load_ps(&in[y_w + x]);
                        // load pixels to the right, using shuffles instead of unaligned access
                        __m128 right = _mm_load_ss(&in[y_w + x + 4]);
                        right = _mm_move_ss(here, right);
                        right = _mm_shuffle_ps(right, right, (0 << 6) + (3 << 4) + (2 << 2) + (1 << 0));
                        // difference
                        __m128 g_x = right - here;

                        // save
                        _mm_store_ps(&in_x[y_w + x], g_x);

                        // load pixels below
                        __m128 down = _mm_load_ps(&in[y_w + w + x]);
                        // difference
                        __m128 g_y = down - here;
                        // save
                        _mm_store_ps(&in_y[y_w + x], g_y);

                        // norm
                        __m128 gnorm = _mm_sqrt_ps(g_x * g_x + g_y * g_y);
                        // add norm to total
                        tv += gnorm[0];
                        tv += gnorm[1];
                        tv += gnorm[2];
                        tv += gnorm[3];

                        // set zeroes to infinity
                        gnorm = _mm_or_ps(gnorm, _mm_and_ps(mm_inf, _mm_cmpeq_ps(gnorm, mm_zero)));

                        // compute gradient, using shuffling to get right gradient working
                        // right
                        __m128 objective_right = _mm_load_ss(&objective_gradient[y_w + x + 4]);
                        __m128 gradient_right = g_x / gnorm;
                        gradient_right = _mm_shuffle_ps(gradient_right, gradient_right, (2 << 6) + (1 << 4) + (0 << 2) + (3 << 0));
                        objective_right = _mm_add_ss(objective_right, gradient_right);
                        _mm_store_ss(&objective_gradient[y_w + x + 4], objective_right);

                        // here
                        gradient_right = _mm_move_ss(gradient_right, mm_zero);
                        __m128 objective_here = _mm_load_ps(&objective_gradient[y_w + x]);
                        objective_here += gradient_right;
                        objective_here += -(g_x + g_y) / gnorm;
                        _mm_store_ps(&objective_gradient[y_w + x], objective_here);

                        // bottom
                        __m128 objective_down = _mm_load_ps(&objective_gradient[y_w + w + x]);
                        _mm_store_ps(&objective_gradient[y_w + w + x], objective_down + g_y / gnorm);
                }
                // last 3 normal pixels
                for(unsigned x = simd_w; x < w-1; x++) {
                        float g_x = *p(in, x+1, y, w, h) - *p(in, x, y, w, h);
                        *p(in_x, x, y, w, h) = g_x;
                        float g_y = *p(in, x, y+1, w, h) - *p(in, x, y, w, h);
                        *p(in_y, x, y, w, h) = g_y;
                        float g_norm = sqrt(sqr(g_x) + sqr(g_y));
                        tv += g_norm;
                        if(g_norm != 0) {
                                *p(objective_gradient, x, y, w, h) += -(g_x + g_y) / g_norm;
                                *p(objective_gradient, x+1, y, w, h) += g_x / g_norm;
                                *p(objective_gradient, x, y+1, w, h) += g_y / g_norm;
                        }
                }
                // edge column pixel
                *p(in_x, w-1, y, w, h) = 0.;
                float g_y = *p(in, w-1, y+1, w, h) - *p(in, w-1, y, w, h);
                *p(in_y, w-1, y, w, h) = g_y;
                float g_norm = sqrt(sqr(g_y));
                tv += g_norm;
                if(g_norm != 0) {
                        *p(objective_gradient, w-1, y, w, h) += -(g_y) / g_norm;
                        *p(objective_gradient, w-1, y+1, w, h) += g_y / g_norm;
                }
        }
        // edge row
        for(unsigned x = 0; x < w-1; x++) {
                float g_x = *p(in, x+1, h-1, w, h) - *p(in, x, h-1, w, h);
                *p(in_x, x, h-1, w, h) = g_x;
                *p(in_y, x, h-1, w, h) = 0.;
                float g_norm = sqrt(sqr(g_x));
                tv += g_norm;
                if(g_norm != 0) {
                        *p(objective_gradient, x, h-1, w, h) += -(g_x) / g_norm;
                        *p(objective_gradient, x+1, h-1, w, h) += g_x / g_norm;
                }
        }
        // corner pixel has zero gradient, don't need to do anything
        *p(in_x, w-1, h-1, w, h) = 0.;
        *p(in_y, w-1, h-1, w, h) = 0.;
        return tv;
}
