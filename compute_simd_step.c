#include <xmmintrin.h>
#include <math.h>
#include "utils.h"

// SSE2, optimized versions of functions in compute.c

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
        const __m128 minf = _mm_set_ps1(INFINITY);
        const __m128 mzero = _mm_set_ps1(0.);

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
        __m128 g_norm = mzero;
        for(unsigned c = 0; c < nchannel; c++) {
                g_norm += SQR(g_xs[c]);
                g_norm += SQR(g_ys[c]);
        }
        g_norm = _mm_sqrt_ps(g_norm);

        float alpha = 1./sqrtf(nchannel);
        *tv += alpha * g_norm[0];
        *tv += alpha * g_norm[1];
        *tv += alpha * g_norm[2];
        *tv += alpha * g_norm[3];

        __m128 malpha = _mm_set_ps1(alpha);

        // set zeroes to infinity
        g_norm = _mm_or_ps(g_norm, _mm_and_ps(minf, _mm_cmpeq_ps(g_norm, mzero)));

        // compute derivatives
        for(unsigned c = 0; c < nchannel; c++) {
                __m128 g_x = g_xs[c];
                __m128 g_y = g_ys[c];
                struct aux *aux = &auxs[c];

                // N.B. for numerical stability and same exact result as the c version,
                // we must calculate the objective gradient at x+1 before x
                {
                        float *pobj_r = p(aux->obj_gradient, x+1, y, w, h);
                        __m128 obj_r = _mm_loadu_ps(pobj_r);
                        obj_r += malpha * g_x / g_norm;
                        _mm_storeu_ps(pobj_r, obj_r);
                }

                {
                        float *pobj = p(aux->obj_gradient, x, y, w, h);
                        __m128 obj = _mm_load_ps(pobj);
                        obj += malpha * -(g_x + g_y) / g_norm;
                        _mm_store_ps(pobj, obj);
                }

                {
                        float *pobj_b = p(aux->obj_gradient, x, y+1, w, h);
                        __m128 obj_b = _mm_load_ps(pobj_b);
                        obj_b += malpha * g_y / g_norm;
                        _mm_store_ps(pobj_b, obj_b);
                }
        }
        // store
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                _mm_store_ps(p(aux->temp[0], x, y, w, h), g_xs[c]);
                _mm_store_ps(p(aux->temp[1], x, y, w, h), g_ys[c]);
        }
}

static double compute_step_tv_simd(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel]) {
        if(w < 4) {
                return compute_step_tv_c(w, h, nchannel, auxs);
        }

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

static void clamp_dct_simd(struct coef *coef, float *boxed, unsigned blocks) {
        __m128 mhalf = _mm_set_ps1(0.5);
        for(unsigned i = 0; i < blocks; i++) {
                for(unsigned j = 0; j < 64; j+=4) {
                        __m128 coef_data = _mm_cvtpi16_ps(*(__m64 *)&(coef->data[i*64+j]));
                        __m128 coef_quant_table = _mm_cvtpi16_ps(*(__m64 *)&(coef->quant_table[j]));

                        __m128 min = (coef_data - mhalf) * coef_quant_table;
                        __m128 max = (coef_data + mhalf) * coef_quant_table;
                        __m128 data = _mm_load_ps(&boxed[i*64+j]);
                        data =_mm_max_ps(min, _mm_min_ps(max, data));
                        _mm_store_ps(&boxed[i*64+j], data);
                }
        }
        _mm_empty();
}

static void compute_step_tv2_inner_simd(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel], float alpha, unsigned x, unsigned y, double *tv2) {
        __m128 g_xxs[3] = {0};
        __m128 g_xy_syms[3] = {0};
        __m128 g_yys[3] = {0};

        const __m128 mtwo = _mm_set_ps1(2.);
        const __m128 minf = _mm_set_ps1(INFINITY);
        const __m128 mzero = _mm_set_ps1(0.);

        __m128 malpha = _mm_set_ps1(alpha * 1./sqrtf(nchannel));

        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];

                __m128 g_x = _mm_load_ps(p(aux->temp[0], x, y, w, h));
                __m128 g_y = _mm_load_ps(p(aux->temp[1], x, y, w, h));

                // backward x
                g_xxs[c] = g_x - _mm_loadu_ps(p(aux->temp[0], x-1, y, w, h));
                // backward x
                __m128 g_yx = g_y - _mm_loadu_ps(p(aux->temp[1], x-1, y, w, h));
                // backward y
                __m128 g_xy = g_x - _mm_load_ps(p(aux->temp[0], x, y-1, w, h));
                // backward y
                g_yys[c] = g_y - _mm_load_ps(p(aux->temp[1], x, y-1, w, h));
                // symmetrize
                g_xy_syms[c] = (g_xy + g_yx) / mtwo;
        }

        // norm
        __m128 g2_norm = mzero;
        for(unsigned c = 0; c < nchannel; c++) {
                g2_norm += SQR(g_xxs[c]) + mtwo * SQR(g_xy_syms[c]) + SQR(g_yys[c]);
        }
        g2_norm = _mm_sqrt_ps(g2_norm);

        __m128 alpha_norm = malpha * g2_norm;
        *tv2 += alpha_norm[0];
        *tv2 += alpha_norm[1];
        *tv2 += alpha_norm[2];
        *tv2 += alpha_norm[3];

        // set zeroes to infinity
        g2_norm = _mm_or_ps(g2_norm, _mm_and_ps(minf, _mm_cmpeq_ps(g2_norm, mzero)));

        for(unsigned c = 0; c < nchannel; c++) {
                __m128 g_xx = g_xxs[c];
                __m128 g_yy = g_yys[c];
                __m128 g_xy_sym = g_xy_syms[c];
                struct aux *aux = &auxs[c];

                // N.B. for same exact result as the c version,
                // we must calculate the objective gradient from right to left
                {
                        float *pobj_ur = p(aux->obj_gradient, x+1, y-1, w, h);
                        __m128 obj_ur = _mm_loadu_ps(pobj_ur);
                        obj_ur += malpha * ((-g_xy_sym) / g2_norm);
                        _mm_storeu_ps(pobj_ur, obj_ur);
                }

                {
                        float *pobj_r = p(aux->obj_gradient, x+1, y, w, h);
                        __m128 obj_r = _mm_loadu_ps(pobj_r);
                        obj_r += malpha * ((g_xy_sym + g_xx) / g2_norm);
                        _mm_storeu_ps(pobj_r, obj_r);
                }

                {
                        float *pobj_u = p(aux->obj_gradient, x, y-1, w, h);
                        __m128 obj_u = _mm_load_ps(pobj_u);
                        obj_u += malpha * ((g_yy + g_xy_sym) / g2_norm);
                        _mm_store_ps(pobj_u, obj_u);
                }

                {
                        float *pobj = p(aux->obj_gradient, x, y, w, h);
                        __m128 obj = _mm_load_ps(pobj);
                        obj += malpha * (-(mtwo * g_xx + mtwo * g_xy_sym + mtwo * g_yy) / g2_norm);
                        _mm_store_ps(pobj, obj);
                }

                {
                        float *pobj_b = p(aux->obj_gradient, x, y+1, w, h);
                        __m128 obj_b = _mm_load_ps(pobj_b);
                        obj_b += malpha * ((g_yy + g_xy_sym) / g2_norm);
                        _mm_store_ps(pobj_b, obj_b);
                }

                {
                        float *pobj_l = p(aux->obj_gradient, x-1, y, w, h);
                        __m128 obj_l = _mm_loadu_ps(pobj_l);
                        obj_l += malpha * ((g_xy_sym + g_xx) / g2_norm);
                        _mm_storeu_ps(pobj_l, obj_l);
                }

                {
                        float *pobj_lb = p(aux->obj_gradient, x-1, y+1, w, h);
                        __m128 obj_lb = _mm_loadu_ps(pobj_lb);
                        obj_lb += malpha * ((-g_xy_sym) / g2_norm);
                        _mm_storeu_ps(pobj_lb, obj_lb);
                }
        }
}

static double compute_step_tv2_simd(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel], float alpha) {
        if(w < 8 || h < 2) {
                return compute_step_tv2_c(w, h, nchannel, auxs, alpha);
        }

        double tv2 = 0.;
        for(unsigned x = 0; x < w; x++) {
                compute_step_tv2_inner_c(w, h, nchannel, auxs, alpha, x, 0, &tv2);
        }
        for(unsigned y = 1; y < h-1; y++) {
                for(unsigned x = 0; x < 4; x++) {
                        compute_step_tv2_inner_c(w, h, nchannel, auxs, alpha, x, y, &tv2);
                }
                for(unsigned x = 4; x < w-4; x+=4) {
                        compute_step_tv2_inner_simd(w, h, nchannel, auxs, alpha, x, y, &tv2);
                }
                for(unsigned x = w-4; x < w; x++) {
                        compute_step_tv2_inner_c(w, h, nchannel, auxs, alpha, x, y, &tv2);
                }
        }
        for(unsigned x = 0; x < w; x++) {
                compute_step_tv2_inner_c(w, h, nchannel, auxs, alpha, x, h-1, &tv2);
        }
        return tv2;
}
