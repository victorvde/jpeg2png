#include <xmmintrin.h>
#include <math.h>
#include "utils.h"

static double compute_step_prob_simd(unsigned w, unsigned h, float alpha, struct coef *coef, float *cos, float *obj_gradient) {
        double prob_dist = 0.;
        unsigned block_w = coef->w / 8;
        unsigned block_h = coef->h / 8;
        for(unsigned block_y = 0; block_y < block_h; block_y++) {
                for(unsigned block_x = 0; block_x < block_w; block_x++) {
                        unsigned i = block_y * block_w + block_x;
                        float *cosb = &cos[i*64];
                        for(unsigned j = 0; j < 64; j+=4) {
                                __m128 coef_data = _mm_cvtpi16_ps(*(__m64 *)&(coef->data[i*64+j]));
                                __m128 coef_quant_table = _mm_cvtpi16_ps(*(__m64 *)&(coef->quant_table[j]));
                                _mm_empty();

                                __m128 cosb_j = _mm_load_ps(&cosb[j]);
                                cosb_j = cosb_j - coef_data * coef_quant_table;
                                __m128 dist = SQR(cosb_j / coef_quant_table);
                                prob_dist += dist[0];
                                prob_dist += dist[1];
                                prob_dist += dist[2];
                                prob_dist += dist[3];
                                cosb_j = cosb_j / SQR(coef_quant_table);
                                _mm_store_ps(&cosb[j], cosb_j);
                        }
                        idct8x8s(cosb);
                        if(coef->w_samp > 1 || coef->h_samp > 1) {
                                 for(unsigned in_y = 0; in_y < 8; in_y++) {
                                        for(unsigned in_x = 0; in_x < 8; in_x++) {
                                                unsigned j = in_y * 8 + in_x;
                                                unsigned cx = block_x * 8 + in_x;
                                                unsigned cy = block_y * 8 + in_y;
                                                for(unsigned sy = 0; sy < coef->h_samp; sy++) {
                                                        for(unsigned sx = 0; sx < coef->w_samp; sx++) {
                                                                unsigned y = cy * coef->h_samp + sy;
                                                                unsigned x = cx * coef->w_samp + sx;
                                                                *p(obj_gradient, x, y, w, h) += alpha * cosb[j];
                                                        }
                                                }
                                        }
                                }
                        } else {
                                __m128 malpha = _mm_set_ps1(alpha);
                                for(unsigned j = 0; j < 64; j+=4) {
                                        unsigned in_y = j / 8;
                                        unsigned in_x = j % 8;
                                        unsigned x = block_x * 8 + in_x;
                                        unsigned y = block_y * 8 + in_y;
                                        __m128 obj = _mm_load_ps(&obj_gradient[y*w+x]);
                                        __m128 cosb_j = _mm_load_ps(&cosb[j]);
                                        obj += malpha * cosb_j;
                                        _mm_store_ps(&obj_gradient[y*w+x], obj);
                                }
                        }
                }
        }
        return 0.5 * prob_dist;
}

static double compute_step_tv_simd(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel]) {
        return compute_step_tv_c(w, h, nchannel, auxs);
}
/*
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
*/
