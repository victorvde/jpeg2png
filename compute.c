#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "jpeg2png.h"
#include "compute.h"
#include "utils.h"
#include "box.h"
#include "logger.h"

#include "ooura/dct.h"

static_assert(FLT_EVAL_METHOD == 0, "to preserve identical output please disable excess precision");
#ifdef PRAGMA_FP_CONTRACT
#pragma STDC FP_CONTRACT OFF
#endif

struct aux {
        float *obj_gradient;
        // temp[0] is (TVG) d_xs (PROJ/PROB) boxed
        // temp[1] is (TVG) d_ys (PROJ/PROB) subsampled
        // temp[2] is (TVG) d_norms (PROJ/PROB) DCTed data
        float *temp[3];
        float *fdata;
        float *fista;
};

// Note: destroys cos
POSSIBLY_UNUSED static double compute_step_prob_c(unsigned w, unsigned h, float alpha, struct coef *coef, float *cos, float *obj_gradient) {
        double prob_dist = 0.;
        unsigned block_w = coef->w / 8;
        unsigned block_h = coef->h / 8;
        for(unsigned block_y = 0; block_y < block_h; block_y++) {
                for(unsigned block_x = 0; block_x < block_w; block_x++) {
                        unsigned i = block_y * block_w + block_x;
                        float *cosb = &cos[i*64];
                        for(unsigned j = 0; j < 64; j++) {
                                cosb[j] -= (float)coef->data[i*64+j] * coef->quant_table[j];
                                prob_dist += 0.5 * sqr(cosb[j] / coef->quant_table[j]);
                                cosb[j] = cosb[j] / sqr((float)coef->quant_table[j]);
                        }
                        idct8x8s(cosb);
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
                }
        }
        return alpha * prob_dist;
}

static void compute_step_tv_norms_c(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel], unsigned x, unsigned y, float *d_norms) {
        float d_norm = 0.;
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                // forward difference x
                float d_x = x >= w-1 ? 0. : *p(aux->fdata, x+1, y, w, h) - *p(aux->fdata, x, y, w, h);
                // forward difference y
                float d_y = y >= h-1 ? 0. : *p(aux->fdata, x, y+1, w, h) - *p(aux->fdata, x, y, w, h);
                // store
                *p(aux->temp[0], x, y, w, h) = d_x;
                *p(aux->temp[1], x, y, w, h) = d_y;
                // norm
                d_norm += sqr(d_x);
                d_norm += sqr(d_y);
        }
        d_norm = sqrt(d_norm);
        *p(d_norms, x, y, w, h) = d_norm;
}

static void compute_step_tv_objective_c(unsigned w, unsigned h, float alpha, struct aux *aux, unsigned x, unsigned y, float *d_norms) {
        // compute derivatives
        float obj = 0.;

        float d_norm = *p(d_norms, x, y, w, h);
        if(d_norm != 0.) {
                // forward difference x
                float d_x = *p(aux->temp[0], x, y, w, h);
                // forward difference y
                float d_y = *p(aux->temp[1], x, y, w, h);
                obj += -(d_x + d_y) / d_norm;
        }

        if(x > 0) {
                float dnorm_r = *p(d_norms, x-1, y, w, h);
                if(dnorm_r != 0.) {
                        // forward difference x (other pixel)
                        float d_x_o = *p(aux->temp[0], x-1, y, w, h);
                        obj += d_x_o / dnorm_r;
                        }
        }

        if(y > 0) {
                float dnorm_u = *p(d_norms, x, y-1, w, h);
                if(dnorm_u != 0.) {
                        // forward difference y (other pixel)
                        float d_y_o = *p(aux->temp[1], x, y-1, w, h);
                        obj += d_y_o / dnorm_u;
                }
        }

        *p(aux->obj_gradient, x, y, w, h) += alpha * obj;
}

static double compute_step_tv_c(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel]) {
        double tv = 0.;
        float alpha = 1./sqrt(nchannel);
        float *d_norms = auxs[0].temp[2];

        OPENMP(parallel for collapse(2))
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        compute_step_tv_norms_c(w, h, nchannel, auxs, x, y, d_norms);
                }
        }

        OPENMP(parallel for collapse(2))
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        float d_norm = *p(d_norms, x, y, w, h);
                        tv += alpha * d_norm;
                }
        }

        OPENMP(parallel for collapse(3))
        for(unsigned c = 0; c < nchannel; c++) {
                for(unsigned y = 0; y < h; y++) {
                        for(unsigned x = 0; x < w; x++) {
                                struct aux *aux = &auxs[c];
                                compute_step_tv_objective_c(w, h, alpha, aux, x, y, d_norms);
                        }
                }
        }
        return tv;
}

static void compute_step_tv2_inner_c(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel], float alpha, unsigned x, unsigned y, double *tv2) {
        float g_xxs[3] = {0};
        float g_xy_syms[3] = {0};
        float g_yys[3] = {0};

        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];

                // backward x
                g_xxs[c] = x <= 0 ? 0. : *p(aux->temp[0], x, y, w, h) - *p(aux->temp[0], x-1, y, w, h);
                // backward x
                float g_yx = x <= 0 ? 0. : *p(aux->temp[1], x, y, w, h) - *p(aux->temp[1], x-1, y, w, h);
                // backward y
                float g_xy = y <= 0 ? 0. : *p(aux->temp[0], x, y, w, h) - *p(aux->temp[0], x, y-1, w, h);
                // backward y
                g_yys[c] = y <= 0 ? 0. : *p(aux->temp[1], x, y, w, h) - *p(aux->temp[1], x, y-1, w, h);
                // symmetrize
                g_xy_syms[c] = (g_xy + g_yx) / 2.;
        }
        // norm
        float g2_norm = 0.;
        for(unsigned c = 0; c < nchannel; c++) {
                g2_norm += sqr(g_xxs[c]) + 2 * sqr(g_xy_syms[c]) + sqr(g_yys[c]);
        }
        g2_norm = sqrt(g2_norm);

        alpha = alpha * 1./sqrt(nchannel);
        *tv2 += alpha * g2_norm;

        // compute derivatives
        if(g2_norm != 0.) {
                for(unsigned c = 0; c < nchannel; c++) {
                        float g_xx = g_xxs[c];
                        float g_yy = g_yys[c];
                        float g_xy_sym = g_xy_syms[c];
                        struct aux *aux = &auxs[c];

                        *p(aux->obj_gradient, x, y, w, h) += alpha * (-(2 * g_xx + 2 * g_xy_sym + 2 * g_yy) / g2_norm);
                        if(x > 0) {
                                *p(aux->obj_gradient, x-1, y, w, h) += alpha * ((g_xy_sym + g_xx) / g2_norm);
                        }
                        if(x < w-1) {
                                *p(aux->obj_gradient, x+1, y, w, h) += alpha * ((g_xy_sym + g_xx) / g2_norm);
                        }
                        if(y > 0) {
                                *p(aux->obj_gradient, x, y-1, w, h) += alpha * ((g_yy + g_xy_sym) / g2_norm);
                        }
                        if(y < h-1) {
                                *p(aux->obj_gradient, x, y+1, w, h) += alpha * ((g_yy + g_xy_sym) / g2_norm);
                        }
                        if(x < w-1 && y > 0) {
                                *p(aux->obj_gradient, x+1, y-1, w, h) += alpha * ((-g_xy_sym) / g2_norm);
                        }
                        if(x > 0 && y < h-1) {
                                *p(aux->obj_gradient, x-1, y+1, w, h) += alpha * ((-g_xy_sym) / g2_norm);
                        }
                }
        }
}

static double compute_step_tv2_c(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel], float alpha) {
        double tv2 = 0.;
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        compute_step_tv2_inner_c(w, h, nchannel, auxs, alpha, x, y, &tv2);
                }
        }
        return tv2;
}

static double compute_norm(unsigned w, unsigned h, float *data) {
        double norm = 0.;
        for(unsigned i = 0; i < h * w; i++) {
                norm += sqr(data[i]);
        }
        return sqrt(norm);
}

static void compute_do_step(unsigned w, unsigned h, float *fdata, float *obj_gradient, float step_size) {
        // do step (normalized)
        float norm = compute_norm(w, h, obj_gradient);
        if(norm != 0.) {
                for(unsigned i = 0; i < h * w; i++) {
                        fdata[i] = fdata[i] - step_size * (obj_gradient[i] /  norm);
                }
        }
}

#ifdef USE_SIMD
#include "compute_simd_step.c"
#endif

static double compute_step(
        unsigned w, unsigned h,
        unsigned nchannel,
        struct coef coefs[nchannel], struct aux auxs[nchannel],
        float step_size, float weight, float pweight[nchannel],
        struct logger *log)
{
        float total_alpha = 0.;

        double prob_dist = 0.;
        OPENMP(parallel for reduction(+:total_alpha) reduction(+:prob_dist))
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                struct coef *coef = &coefs[c];

                // initialize gradient
                for(unsigned i = 0; i < h * w; i++) {
                        aux->obj_gradient[i] = 0.;
                }

                // DCT coefficent distance
                if(pweight[c] !=  0.) {
                        float p_alpha = pweight[c] * 2 * 255 * sqrt(2);
                        total_alpha += p_alpha;
                        prob_dist += POSSIBLY_SIMD(compute_step_prob)(w, h, p_alpha, coef, aux->temp[2], aux->obj_gradient);
                }
        }

        // TV
        total_alpha += nchannel;
        double tv = POSSIBLY_SIMD(compute_step_tv)(w, h, nchannel, auxs);

        // TVG second order
        double tv2 = 0.;
        if(weight != 0.) {
                float alpha = weight / sqrt(4 / 2);
                total_alpha += alpha * nchannel;
                tv2 = POSSIBLY_SIMD(compute_step_tv2)(w, h, nchannel, auxs, alpha);
        }

        // do step
        OPENMP(parallel for)
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                compute_do_step(w, h, aux->fdata, aux->obj_gradient, step_size);
        }

        double objective = (tv + tv2 + prob_dist) / total_alpha;
        logger_log(log, objective, prob_dist, tv, tv2);

        return objective;
}

static void aux_init(unsigned w, unsigned h, struct coef *coef, struct aux *aux) {
        for(unsigned i = 0; i < 3; i++) {
                float *t = alloc_real(h * w);
                aux->temp[i] = t;
        }
        float *cos = aux->temp[2];
        unsigned blocks = (coef->h / 8) * (coef->w / 8);
        for(unsigned i = 0; i < blocks; i++) {
                for(unsigned j = 0; j < 64; j++) {
                        cos[i*64+j] = coef->data[i*64+j] * coef->quant_table[j];
                }
        }

        float *obj_gradient = alloc_real(h * w);
        aux->obj_gradient = obj_gradient;

        float *fdata = alloc_real(h * w);
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        unsigned cy = MIN(y / coef->h_samp, coef->h-1);
                        unsigned cx = MIN(x / coef->w_samp, coef->w-1);
                        *p(fdata, x, y, w, h) = *p(coef->fdata, cx, cy, coef->w, coef->h);
                }
        }
        aux->fdata = fdata;

        float *fista = alloc_real(h * w);
        memcpy(fista, fdata, sizeof(float) * w * h);
        aux->fista = fista;
}

static void aux_destroy(struct aux *aux) {
        for(unsigned i = 0; i < 3; i++) {
                free_real(aux->temp[i]);
        }
        free_real(aux->obj_gradient);
        free_real(aux->fdata);
        free_real(aux->fista);
}

POSSIBLY_UNUSED static void clamp_dct_c(struct coef *coef, float *boxed, unsigned blocks) {
                for(unsigned i = 0; i < blocks; i++) {
                for(unsigned j = 0; j < 64; j++) {
                        float min = (coef->data[i*64+j] - 0.5f) * coef->quant_table[j];
                        float max = (coef->data[i*64+j] + 0.5f) * coef->quant_table[j];
                        boxed[i*64+j] = CLAMP(boxed[i*64+j], min, max);
                }
        }
}

static void compute_projection(unsigned w, unsigned h, struct aux *aux, struct coef *coef) {
        unsigned blocks = (coef->h / 8) * (coef->w / 8);
        float *subsampled;
        float *boxed = aux->temp[0];
        bool resample = !(coef->w == w && coef->h == h);

        if(resample) {
                subsampled = aux->temp[1];
        } else {
                subsampled = aux->fdata;
        }

        if(resample) {
                for(unsigned cy = 0; cy < coef->h; cy++) {
                        for(unsigned cx = 0; cx < coef->w; cx++) {
                                float mean = 0.;
                                for(unsigned sy = 0; sy < coef->h_samp; sy++) {
                                        for(unsigned sx = 0; sx < coef->w_samp; sx++) {
                                                unsigned y = cy * coef->h_samp + sy;
                                                unsigned x = cx * coef->w_samp + sx;
                                                mean += *p(aux->fdata, x, y, w, h);
                                        }
                                }
                                mean /= coef->w_samp * coef->h_samp;
                                *p(subsampled, cx, cy, coef->w, coef->h) = mean;
                                for(unsigned sy = 0; sy < coef->h_samp; sy++) {
                                        for(unsigned sx = 0; sx < coef->w_samp; sx++) {
                                                unsigned y = cy * coef->h_samp + sy;
                                                unsigned x = cx * coef->w_samp + sx;
                                                *p(aux->fdata, x, y, w, h) -= mean;
                                }
                                }
                        }
                }
        }

        box(subsampled, boxed, coef->w, coef->h);

        for(unsigned i = 0; i < blocks; i++) {
                dct8x8s(&boxed[i*64]);
        }

        POSSIBLY_SIMD(clamp_dct)(coef, boxed, blocks);

        memcpy(aux->temp[2], boxed, coef->w * coef->h * sizeof(float));

        for(unsigned i = 0; i < blocks; i++) {
                idct8x8s(&boxed[i*64]);
        }

        unbox(boxed, subsampled, coef->w, coef->h);

        if(resample) {
                for(unsigned cy = 0; cy < coef->h; cy++) {
                        for(unsigned cx = 0; cx < coef->w; cx++) {
                                float mean = *p(subsampled, cx, cy, coef->w, coef->h);
                                for(unsigned sy = 0; sy < coef->h_samp; sy++) {
                                        for(unsigned sx = 0; sx < coef->w_samp; sx++) {
                                                unsigned y = cy * coef->h_samp + sy;
                                                unsigned x = cx * coef->w_samp + sx;
                                                *p(aux->fdata, x, y, w, h) += mean;
                                        }
                                }
                        }
                }
        }
}

void compute(unsigned nchannel, struct coef coefs[nchannel], struct logger *log, struct progressbar *pb, float weight, float pweight[nchannel], unsigned iterations) {
        unsigned h = 0;
        unsigned w = 0;
        for(unsigned c = 0; c < nchannel; c++) {
                struct coef *coef = &coefs[c];
                w = MAX(w, coef->w * coef->w_samp);
                h = MAX(h, coef->h * coef->h_samp);
        }
        ASSUME(w % 8 == 0);
        ASSUME(h % 8 == 0);
        struct aux *auxs = malloc(sizeof(*auxs) * nchannel);
        for(unsigned c = 0; c < nchannel; c++) {
                aux_init(w, h, &coefs[c], &auxs[c]);
        }

        float radius = sqrt(w*h) / 2;
        float t = 1;
        for(unsigned i = 0; i < iterations; i++) {
                log->iteration = i;

                float tnext = (1 + sqrt(1 + 4 * sqr(t))) / 2;
                float factor = (t - 1) / tnext;
                for(unsigned c = 0; c < nchannel; c++) {
                        struct aux *aux = &auxs[c];
                        for(unsigned j = 0; j < w * h; j++) {
                                aux->fista[j] = aux->fdata[j] + factor * (aux->fdata[j] - aux->fista[j]);
                        }
                        SWAP(float *, aux->fdata, aux->fista);
                }
                t = tnext;

                compute_step(w, h, nchannel, coefs, auxs, radius / sqrt(1 + iterations), weight, pweight, log);
                OPENMP(parallel for schedule(dynamic))
                for(unsigned c = 0; c < nchannel; c++) {
                        compute_projection(w, h, &auxs[c], &coefs[c]);
                }
                if(pb) {
                        OPENMP(critical(progressbar))
                        progressbar_inc(pb);
                }
        }

        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                struct coef *coef = &coefs[c];
                SWAP(float *, aux->fdata, coef->fdata);
                coef->w = w;
                coef->h = h;
                aux_destroy(aux);
        }
        free(auxs);
}
