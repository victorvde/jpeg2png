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

// working buffers for each component
struct aux {
        // DCT coefficients for step_prob
        float *cos;
        // gradient (derivative) of the objective function
        float *obj_gradient;
        // temp[0] = pixel differences in x direction
        // temp[1] = pixel differences in y direction
        // also used differently in compute_projection
        float *temp[2];
        // image data
        float *fdata;
        // previous step image data for FISTA
        float *fista;
};

// compute objective gradient for the distance of DCT coefficients from normal decoding
// N.B. destroys cos
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
                                prob_dist += 0.5 * sqr(cosb[j] / coef->quant_table[j]); // objective function
                                cosb[j] = cosb[j] / sqr((float)coef->quant_table[j]); // derivative
                        }
                        idct8x8s(cosb);
                        // unbox and possibly upsample derivative
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

// compute objective gradient for TV for one pixel
static void compute_step_tv_inner_c(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel], unsigned x, unsigned y, double *tv) {
        float g_xs[3] = {0};
        float g_ys[3] = {0};
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                // forward difference x
                g_xs[c] = x >= w-1 ? 0. : *p(aux->fdata, x+1, y, w, h) - *p(aux->fdata, x, y, w, h);
                // forward difference y
                g_ys[c] = y >= h-1 ? 0. : *p(aux->fdata, x, y+1, w, h) - *p(aux->fdata, x, y, w, h);
        }
        // norm
        float g_norm = 0.;
        for(unsigned c = 0; c < nchannel; c++) {
                g_norm += sqr(g_xs[c]);
                g_norm += sqr(g_ys[c]);
        }
        g_norm = sqrtf(g_norm);
        float alpha = 1./sqrtf(nchannel);
        *tv += alpha * g_norm; // objective function
        // compute derivatives (see notes)
        for(unsigned c = 0; c < nchannel; c++) {
                float g_x = g_xs[c];
                float g_y = g_ys[c];
                struct aux *aux = &auxs[c];
                if(g_norm != 0) {
                        *p(aux->obj_gradient, x, y, w, h) += alpha * -(g_x + g_y) / g_norm;
                        if(x < w-1) {
                                *p(aux->obj_gradient, x+1, y, w, h) += alpha * g_x / g_norm;
                        }
                        if(y < h-1) {
                                *p(aux->obj_gradient, x, y+1, w, h) += alpha * g_y / g_norm;
                        }
                }
        }
        // store for use in tv2
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                *p(aux->temp[0], x, y, w, h) = g_xs[c];
                *p(aux->temp[1], x, y, w, h) = g_ys[c];
        }
}

// compute objective gradient for TV
static double compute_step_tv_c(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel]) {
        double tv = 0.;
        ASSUME(nchannel <= 3);
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        compute_step_tv_inner_c(w, h, nchannel, auxs, x, y, &tv);
                }
        }
        return tv;
}

// compute objective gradient for second order TGV for one pixel
static void compute_step_tv2_inner_c(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel], float alpha, unsigned x, unsigned y, double *tv2) {
        float g_xxs[3] = {0};
        float g_xy_syms[3] = {0};
        float g_yys[3] = {0};

        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];

                // backward difference x
                g_xxs[c] = x <= 0 ? 0. : *p(aux->temp[0], x, y, w, h) - *p(aux->temp[0], x-1, y, w, h);
                // backward difference x
                float g_yx = x <= 0 ? 0. : *p(aux->temp[1], x, y, w, h) - *p(aux->temp[1], x-1, y, w, h);
                // backward difference y
                float g_xy = y <= 0 ? 0. : *p(aux->temp[0], x, y, w, h) - *p(aux->temp[0], x, y-1, w, h);
                // backward difference y
                g_yys[c] = y <= 0 ? 0. : *p(aux->temp[1], x, y, w, h) - *p(aux->temp[1], x, y-1, w, h);
                // symmetrize
                g_xy_syms[c] = (g_xy + g_yx) / 2.;
        }
        // norm
        float g2_norm = 0.;
        for(unsigned c = 0; c < nchannel; c++) {
                g2_norm += sqr(g_xxs[c]) + 2 * sqr(g_xy_syms[c]) + sqr(g_yys[c]);
        }
        g2_norm = sqrtf(g2_norm);

        alpha = alpha * 1./sqrtf(nchannel);
        *tv2 += alpha * g2_norm; // objective function

        // compute derivatives (see notes)
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

// compute objective gradient for second order TGV
static double compute_step_tv2_c(unsigned w, unsigned h, unsigned nchannel, struct aux auxs[nchannel], float alpha) {
        double tv2 = 0.;
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        compute_step_tv2_inner_c(w, h, nchannel, auxs, alpha, x, y, &tv2);
                }
        }
        return tv2;
}

// compute Euclidean norm
static double compute_norm(unsigned w, unsigned h, float *data) {
        double norm = 0.;
        for(unsigned i = 0; i < h * w; i++) {
                norm += sqr(data[i]);
        }
        return sqrtf(norm);
}

// make step in the direction of the objective gradient with distance step_size
static void compute_do_step(unsigned w, unsigned h, float *fdata, float *obj_gradient, float step_size) {
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

// compute objective gradient and make step
static double compute_step(
        unsigned w, unsigned h,
        unsigned nchannel,
        struct coef coefs[nchannel], struct aux auxs[nchannel],
        float step_size, float weight, float pweight[nchannel],
        struct logger *log)
{
        float total_alpha = 0.;

        double prob_dist = 0.;
        OPENMP(parallel for schedule(dynamic) reduction(+:total_alpha) reduction(+:prob_dist))
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                struct coef *coef = &coefs[c];

                // initialize gradient
                for(unsigned i = 0; i < h * w; i++) {
                        aux->obj_gradient[i] = 0.;
                }

                // DCT coefficent distance
                if(pweight[c] !=  0.) {
                        float p_alpha = pweight[c] * 2 * 255 * sqrtf(2);
                        total_alpha += p_alpha;
                        prob_dist += POSSIBLY_SIMD(compute_step_prob)(w, h, p_alpha, coef, aux->cos, aux->obj_gradient);
                }
        }

        // TV
        total_alpha += nchannel;
        double tv = POSSIBLY_SIMD(compute_step_tv)(w, h, nchannel, auxs);

        // TGV second order
        double tv2 = 0.;
        if(weight != 0.) {
                float alpha = weight / sqrtf(4 / 2);
                total_alpha += alpha * nchannel;
                tv2 = POSSIBLY_SIMD(compute_step_tv2)(w, h, nchannel, auxs, alpha);
        }

        // do step
        OPENMP(parallel for schedule(dynamic))
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                compute_do_step(w, h, aux->fdata, aux->obj_gradient, step_size);
        }

        // log objective values
        double objective = (tv + tv2 + prob_dist) / total_alpha;
        logger_log(log, objective, prob_dist, tv, tv2);

        return objective;
}

// initialize working buffers
static void aux_init(unsigned w, unsigned h, struct coef *coef, struct aux *aux) {
        float *cos = alloc_real(coef->h * coef->w);
        unsigned blocks = (coef->h / 8) * (coef->w / 8);
        for(unsigned i = 0; i < blocks; i++) {
                for(unsigned j = 0; j < 64; j++) {
                        cos[i*64+j] = coef->data[i*64+j] * coef->quant_table[j];
                }
        }
        aux->cos = cos;

        for(unsigned i = 0; i < 2; i++) {
                float *t = alloc_real(h * w);
                aux->temp[i] = t;
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
        free_real(coef->fdata);
        coef->fdata = NULL;

        float *fista = alloc_real(h * w);
        memcpy(fista, fdata, sizeof(float) * w * h);
        aux->fista = fista;
}

// destroy working buffers, except the fdata that is returned
static void aux_destroy(struct aux *aux) {
        free_real(aux->cos);
        for(unsigned i = 0; i < 2; i++) {
                free_real(aux->temp[i]);
        }
        free_real(aux->obj_gradient);
        free_real(aux->fista);
}

// clamp the DCT values to interval that quantizes to our jpg
POSSIBLY_UNUSED static void clamp_dct_c(struct coef *coef, float *boxed, unsigned blocks) {
        for(unsigned i = 0; i < blocks; i++) {
                for(unsigned j = 0; j < 64; j++) {
                        float min = (coef->data[i*64+j] - 0.5f) * coef->quant_table[j];
                        float max = (coef->data[i*64+j] + 0.5f) * coef->quant_table[j];
                        boxed[i*64+j] = CLAMP(boxed[i*64+j], min, max);
                }
        }
}

// compute projection of data onto the feasible set defined by our jpg
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

        // downsample and keep the difference
        // more formally, decompose each subsampling block in the direction of our subsampling vector (a vector of ones)
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

        // project onto our DCT box
        box(subsampled, boxed, coef->w, coef->h);

        for(unsigned i = 0; i < blocks; i++) {
                dct8x8s(&boxed[i*64]);
        }

        POSSIBLY_SIMD(clamp_dct)(coef, boxed, blocks);

        memcpy(aux->cos, boxed, coef->w * coef->h * sizeof(float)); // save a copy of the DCT values for step_prob

        for(unsigned i = 0; i < blocks; i++) {
                idct8x8s(&boxed[i*64]);
        }

        unbox(boxed, subsampled, coef->w, coef->h);

        // add back the difference (orthogonal to our subsampling vector)
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

// subgradient method with iteration steps
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
        // working buffers per channel
        struct aux *auxs = malloc(sizeof(*auxs) * nchannel);
        for(unsigned c = 0; c < nchannel; c++) {
                aux_init(w, h, &coefs[c], &auxs[c]);
        }

        float radius = sqrtf(w*h) / 2; // radius of [-0.5, 0.5]^(w*h)
        float t = 1;
        for(unsigned i = 0; i < iterations; i++) {
                log->iteration = i;

                // FISTA
                float tnext = (1 + sqrtf(1 + 4 * sqr(t))) / 2;
                float factor = (t - 1) / tnext;
                for(unsigned c = 0; c < nchannel; c++) {
                        struct aux *aux = &auxs[c];
                        for(unsigned j = 0; j < w * h; j++) {
                                aux->fista[j] = aux->fdata[j] + factor * (aux->fdata[j] - aux->fista[j]);
                        }
                        SWAP(float *, aux->fdata, aux->fista);
                }
                t = tnext;

                // take a step
                compute_step(w, h, nchannel, coefs, auxs, radius / sqrtf(1 + iterations), weight, pweight, log);
                // project back onto feasible set
                OPENMP(parallel for schedule(dynamic))
                for(unsigned c = 0; c < nchannel; c++) {
                        compute_projection(w, h, &auxs[c], &coefs[c]);
                }
                if(pb) {
                        OPENMP(critical(progressbar))
                        progressbar_inc(pb);
                }
        }
        // return result
        for(unsigned c = 0; c < nchannel; c++) {
                struct aux *aux = &auxs[c];
                struct coef *coef = &coefs[c];
                coef->fdata = aux->fdata;
                aux->fdata = NULL;
                coef->w = w;
                coef->h = h;
                aux_destroy(aux);
        }
        free(auxs);
}
