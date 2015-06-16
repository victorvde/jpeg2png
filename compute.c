#include <stdint.h>
#include <string.h>

#include "jpeg2png.h"
#include "compute.h"
#include "utils.h"
#include "box.h"
#include "logger.h"

#include "ooura/dct.h"

struct compute_aux {
        float *q_min;
        float *q_max;
        float *cos;
        float *obj_gradient;
        float *temp[2];
        float *fdata;
        float *fista;
};

// Note: destroys cos
static double compute_step_prob(unsigned w, unsigned h, int16_t *data, float alpha, uint16_t quant_table[64], float *cos, float *obj_gradient) {
        double prob_dist = 0.;
        unsigned block_w = w / 8;
        unsigned block_h = h / 8;
        for(unsigned block_y = 0; block_y < block_h; block_y++) {
                for(unsigned block_x = 0; block_x < block_w; block_x++) {
                        unsigned i = block_y * block_w + block_x;
                        float *cosb = &cos[i*64];
                        for(unsigned j = 0; j < 64; j++) {
                                cosb[j] -= (float)data[i*64+j] * quant_table[j];
                                prob_dist += 0.5 * sqr(cosb[j] / quant_table[j]);
                                cosb[j] = cosb[j] / sqr((float)quant_table[j]);
                        }
                        idct8x8s(cosb);
                        for(unsigned in_y = 0; in_y < 8; in_y++) {
                                for(unsigned in_x = 0; in_x < 8; in_x++) {
                                        unsigned j = in_y * 8 + in_x;
                                        unsigned x = block_x * 8 + in_x;
                                        unsigned y = block_y * 8 + in_y;
                                        *p(obj_gradient, x, y, w, h) += alpha * cosb[j];
                                }
                        }
                }
        }
        return prob_dist;
}

#ifdef USE_SIMD
  #include "compute_simd_step.c"
  #define compute_step_tv compute_step_tv_simd
#else
  #define compute_step_tv compute_step_tv_c
#endif

POSSIBLY_UNUSED static double compute_step_tv_c(unsigned w, unsigned h, float *in, float *obj_gradient, float *in_x, float *in_y) {
        double tv = 0.;
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        // forward gradient x
                        float g_x = x >= w-1 ? 0. : *p(in, x+1, y, w, h) - *p(in, x, y, w, h);
                        // forward gradient y
                        float g_y = y >= h-1 ? 0. : *p(in, x, y+1, w, h) - *p(in, x, y, w, h);
                        // norm
                        float g_norm = sqrt(sqr(g_x) + sqr(g_y));
                        tv += g_norm;
                        // compute derivatives
                        if(g_norm != 0) {
                                *p(obj_gradient, x, y, w, h) += -(g_x + g_y) / g_norm;
                                if(x < w-1) {
                                        *p(obj_gradient, x+1, y, w, h) += g_x / g_norm;
                                }
                                if(y < h-1) {
                                        *p(obj_gradient, x, y+1, w, h) += g_y / g_norm;
                                }
                        }
                        *p(in_x, x, y, w, h) = g_x;
                        *p(in_y, x, y, w, h) = g_y;
                }
        }
        return tv;
}

static double compute_step_tv2(unsigned w, unsigned h, float *obj_gradient, float *in_x, float *in_y, float alpha) {
        double tv2 = 0.;
        for(unsigned y = 0; y < h; y++) {
                for(unsigned x = 0; x < w; x++) {
                        // backward x
                        float g_xx = x <= 0 ? 0. : *p(in_x, x, y, w, h) - *p(in_x, x-1, y, w, h);
                        // backward x
                        float g_yx = x <= 0 ? 0. : *p(in_y, x, y, w, h) - *p(in_y, x-1, y, w, h);
                        // backward y
                        float g_xy = y <= 0 ? 0. : *p(in_x, x, y, w, h) - *p(in_x, x, y-1, w, h);
                        // backward y
                        float g_yy = y <= 0 ? 0. : *p(in_y, x, y, w, h) - *p(in_y, x, y-1, w, h);
                        // symmetrize
                        float g_xy_sym = (g_xy + g_yx) / 2.;
                        // norm
                        float g2_norm = sqrt(sqr(g_xx) + 2 * sqr(g_xy_sym) + sqr(g_yy));
                        tv2 += g2_norm;
                        // compute derivatives
                        if(g2_norm != 0.) {
                                *p(obj_gradient, x, y, w, h) += alpha * (-(2. * g_xx + 2. * g_xy_sym + 2. *  g_yy) / g2_norm);
                                if(x > 0) {
                                        *p(obj_gradient, x-1, y, w, h) += alpha * ((g_xy_sym + g_xx) / g2_norm);
                                }
                                if(x < w-1) {
                                        *p(obj_gradient, x+1, y, w, h) += alpha * ((g_xy_sym + g_xy) / g2_norm);
                                }
                                if(y > 0) {
                                        *p(obj_gradient, x, y-1, w, h) += alpha * ((g_yy + g_xy_sym) / g2_norm);
                                }
                                if(y < h-1) {
                                        *p(obj_gradient, x, y+1, w, h) += alpha * ((g_yy + g_xy_sym) / g2_norm);
                                }
                                if(x < w-1 && y > 0) {
                                        *p(obj_gradient, x+1, y-1, w, h) += alpha * ((-g_xy_sym) / g2_norm);
                                }
                                if(x > 0 && y < h-1) {
                                        *p(obj_gradient, x-1, y+1, w, h) += alpha * ((-g_xy_sym) / g2_norm);
                                }
                        }
                }
        }
        return tv2;
}

static double compute_step(
        unsigned w, unsigned h,
        unsigned ncoef,
        struct coef coefs[ncoef], struct compute_aux auxs[ncoef],
        float step_size, float weight[ncoef], float pweight[ncoef],
        struct logger *log)
{
        float total_alpha = 0.;

        for(unsigned c = 0; c < ncoef; c++) {
                struct compute_aux *aux = &auxs[c];
                for(unsigned i = 0; i < h * w; i++) {
                        aux->obj_gradient[i] = 0.;
                }
        }

        double prob_dist = 0.;
        for(unsigned c = 0; c < ncoef; c++) {
                if(pweight[c] !=  0.) {
                        float p_alpha = pweight[c] * 2. * 255. * sqrt(2.);
                        total_alpha += p_alpha;
                        struct compute_aux *aux = &auxs[c];
                        struct coef *coef = &coefs[c];
                        prob_dist += p_alpha * compute_step_prob(w, h, coef->data, p_alpha, coef->quant_table, aux->cos, aux->obj_gradient);
                }
        }

        double tv = 0.;
        for(unsigned c = 0; c < ncoef; c++) {
                total_alpha += 1.;
                struct compute_aux *aux = &auxs[c];
                tv += compute_step_tv(w, h, aux->fdata, aux->obj_gradient, aux->temp[0], aux->temp[1]);
        }

        double tv2 = 0.;
        for(unsigned c = 0; c < ncoef; c++) {
                if(weight[c] != 0) {
                        float alpha = weight[c] / sqrt(4. / 2.);
                        total_alpha += alpha;
                        struct compute_aux *aux = &auxs[c];
                        tv2 += alpha * compute_step_tv2(w, h, aux->obj_gradient, aux->temp[0], aux->temp[1], alpha);
                }
        }

        for(unsigned c = 0; c < ncoef; c++) {
                struct compute_aux *aux = &auxs[c];

                float norm = 0.;
                for(unsigned i = 0; i < h * w; i++) {
                        norm += sqr(aux->obj_gradient[i]);
                }
                norm = sqrt(norm);

                for(unsigned i = 0; i < h * w; i++) {
                        aux->fdata[i] = aux->fdata[i] - step_size * (aux->obj_gradient[i] /  norm);
                }
        }

        double objective = (tv + tv2 + prob_dist) / total_alpha;
        logger_log(log, objective, prob_dist, tv, tv2);

        return objective;
}

static void compute_aux_init(unsigned w, unsigned h, int16_t *data, uint16_t quant_table[64], float *fdata, struct compute_aux *aux) {
        float *q_max = alloc_real(h * w);
        if(!q_max) { die("allocation error"); }
        float *q_min = alloc_real(h * w);
        if(!q_min) { die("allocation error"); }
        unsigned blocks = (h / 8) * (w / 8);
        for(unsigned i = 0; i < blocks; i++) {
                for(unsigned j = 0; j < 64; j++) {
                        q_max[i*64+j] = (data[i*64+j] + 0.5) * quant_table[j];
                        q_min[i*64+j] = (data[i*64+j] - 0.5) * quant_table[j];
                }
        }
        aux->q_min = q_min;
        aux->q_max = q_max;

        float *cos = alloc_real(h * w);
        if(!cos) { die("allocation error"); }
        for(unsigned i = 0; i < blocks; i++) {
                for(unsigned j = 0; j < 64; j++) {
                        cos[i*64+j] = data[i*64+j] * quant_table[j];
                }
        }
        aux->cos = cos;

        for(unsigned i = 0; i < 2; i++) {
                float *t = alloc_real(h * w);
                if(!t) { die("allocation error"); }
                aux->temp[i] = t;
        }
        float *obj_gradient = alloc_real(h * w);
        if(!obj_gradient) { die("allocation error"); }
        aux->obj_gradient = obj_gradient;

        aux->fdata = fdata;

        float *fista = alloc_real(h * w);
        if(!fista) { die("allocation error"); }
        memcpy(fista, fdata, sizeof(float) * w * h);
        aux->fista = fista;
}

static void compute_aux_destroy(struct compute_aux *aux) {
        free_real(aux->cos);
        free_real(aux->q_min);
        free_real(aux->q_max);
        for(unsigned i = 0; i < 2; i++) {
                free_real(aux->temp[i]);
        }
        free_real(aux->obj_gradient);
        free_real(aux->fista);
}

static void compute_projection(unsigned w, unsigned h, float *fdata, float *temp, float *cos, float *q_min, float *q_max) {
        unsigned blocks = (h / 8) * (w / 8);

        box(fdata, temp, w, h);

        for(unsigned i = 0; i < blocks; i++) {
                dct8x8s(&temp[i*64]);
        }

        for(unsigned i = 0; i < h * w; i++) {
                temp[i] = CLAMP(temp[i], q_min[i], q_max[i]);
        }

        memcpy(cos, temp, w * h * sizeof(float));

        for(unsigned i = 0; i < blocks; i++) {
                idct8x8s(&temp[i*64]);
        }

        unbox(temp, fdata, w, h);
}

void compute(unsigned ncoef, struct coef coefs[ncoef], struct logger *log, struct progressbar *pb, float weight[ncoef], float pweight[ncoef], unsigned iterations) {
        unsigned h = coefs[0].h;
        unsigned w = coefs[0].w;
        for(unsigned c = 1; c < ncoef; c++) {
                ASSUME(coefs[c].w == w);
                ASSUME(coefs[c].h == h);
        }

        struct compute_aux *auxs = malloc(sizeof(*auxs) * ncoef);
        for(unsigned c = 0; c < ncoef; c++) {
                compute_aux_init(w, h, coefs[c].data, coefs[c].quant_table, coefs[c].fdata, &auxs[c]);
        }

        // float radius = sqrt(sqr(w/2.) + sqr(h/2.));
        float radius = sqrt(w*h) / 2;
        float t = 1;
        for(unsigned i = 0; i < iterations; i++) {
                log->iteration = i;

                float tnext = (1 + sqrt(1 + 4 * sqr(t))) / 2;
                float factor = (t - 1) / tnext;
                for(unsigned c = 0; c < ncoef; c++) {
                        struct compute_aux *aux = &auxs[c];
                        for(unsigned j = 0; j < w * h; j++) {
                                aux->fista[j] = aux->fdata[j] + factor * (aux->fdata[j] - aux->fista[j]);
                        }
                        float *f = aux->fdata;
                        aux->fdata = aux->fista;
                        aux->fista = f;
                }
                t = tnext;

                compute_step(w, h, ncoef, coefs, auxs, radius / sqrt(1 + iterations), weight, pweight, log);
                for(unsigned c = 0; c < ncoef; c++) {
                        struct compute_aux *aux = &auxs[c];
                        compute_projection(w, h, aux->fdata, aux->temp[0], aux->cos, aux->q_min, aux->q_max);
                }
                if(pb) {
#ifdef USE_OPENMP
    #pragma omp critical(progressbar)
#endif
                        progressbar_inc(pb);
                }
        }

        for(unsigned c = 0; c < ncoef; c++) {
                struct compute_aux *aux = &auxs[c];
                struct coef *coef = &coefs[c];
                coef->fdata = aux->fdata;
                compute_aux_destroy(aux);
        }
        free(auxs);
}
