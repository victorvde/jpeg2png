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

static void compute_step_tv_inner_simd(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel], unsigned x, unsigned y, double *tv) {
        const __m128 mm_inf = _mm_set_ps1(INFINITY);
        const __m128 mm_zero = _mm_set_ps1(0.);

        __m128 g_xs[3] = {0};
        __m128 g_ys[3] = {0};
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                __m128 here = _mm_load_ps(p(aux->fdata, x, y, w, h));
                // forward gradient x
                g_xs[c] = _mm_loadu_ps(p(aux->fdata, x+1, y, w, h)) - here;
                // forward gradient y
                g_ys[c] = _mm_loadu_ps(p(aux->fdata, x, y+1, w, h)) - here;
        }
        // norm
        __m128 g_norm = mm_zero;
        for(unsigned c = 0; c < nchannel; c++) {
                g_norm += SQR(g_xs[c]);
                g_norm += SQR(g_ys[c]);
        }
        g_norm = _mm_sqrt_ps(g_norm);

        float alpha = 1./sqrt(nchannel);
        *tv += alpha * g_norm[0];
        *tv += alpha * g_norm[1];
        *tv += alpha * g_norm[2];
        *tv += alpha * g_norm[3];

        __m128 malpha = _mm_set_ps1(alpha);

        // set zeroes to infinity
        g_norm = _mm_or_ps(g_norm, _mm_and_ps(mm_inf, _mm_cmpeq_ps(g_norm, mm_zero)));

        // compute derivatives
        for(unsigned c = 0; c < nchannel; c++) {
                __m128 g_x = g_xs[c];
                __m128 g_y = g_ys[c];
                struct aux *aux = &auxs[c];

                // NB: for numerical stability and same exact result as the c version,
                // we must calculate the objective gradient at x+1 before x
                float *pobj_r = p(aux->obj_gradient, x+1, y, w, h);
                __m128 obj_r = _mm_loadu_ps(pobj_r);
                obj_r += malpha * g_x / g_norm;
                _mm_storeu_ps(pobj_r, obj_r);

                float *pobj = p(aux->obj_gradient, x, y, w, h);
                __m128 obj = _mm_load_ps(pobj);
                obj += malpha * -(g_x + g_y) / g_norm;
                _mm_store_ps(pobj, obj);

                float *pobj_b = p(aux->obj_gradient, x, y+1, w, h);
                __m128 obj_b = _mm_load_ps(pobj_b);
                obj_b += malpha * g_y / g_norm;
                _mm_store_ps(pobj_b, obj_b);
        }
        // store
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                _mm_store_ps(p(aux->temp[0], x, y, w, h), g_xs[c]);
                _mm_store_ps(p(aux->temp[1], x, y, w, h), g_ys[c]);
        }
}

static double compute_step_tv_simd(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel]) {
        double tv = 0.;
        ASSUME(nchannel <= 3);
        for(unsigned y = 0; y < h-1; y++) {
                for(unsigned x = 0; x < w-4; x+=4) {
                        compute_step_tv_inner_simd(w, h, nchannel, auxs, x, y, &tv);
                }
                for(unsigned x = w-4; x < w; x++) {
                        compute_step_tv_inner_c(w, h, nchannel, auxs, x, y, &tv);
                }
        }
        for(unsigned x = 0; x < w; x++) {
                compute_step_tv_inner_c(w, h, nchannel, auxs, x, h-1, &tv);
        }
        return tv;
}
