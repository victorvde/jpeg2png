#include <stdint.h>
#include <string.h>

#include "jpeg2png.h"
#include "compute.h"
#include "utils.h"
#include "box.h"
#include "logger.h"

#include "ooura/dct.h"

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

POSSIBLY_UNUSED static void verify_compute_step_tv(unsigned w, unsigned h, double tv, float *in, float *obj_gradient, float *in_x, float *in_y) {
        puts("verify");
        float *obj_gradient_ = alloc_real(w*h);
        float *in_x_ = alloc_real(w*h);
        float *in_y_ = alloc_real(w*h);
        double tv_c = compute_step_tv_c(w, h, in, obj_gradient_, in_x_, in_y_);
        compare("in_x", w, h, in_x_, in_x);
        compare("in_y", w, h, in_y_, in_y);
        compare("obj_gradient", w, h, obj_gradient_, obj_gradient);
        printf("simd %f, original %f\n", tv, tv_c);
        free_real(obj_gradient_);
        free_real(in_y_);
        free_real(in_x_);
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
                        // norm
                        float g2_norm = sqrt(sqr(g_xx) + sqr(g_yx) + sqr(g_xy) + sqr(g_yy));
                        tv2 += g2_norm;
                        // compute derivatives
                        if(g2_norm != 0.) {
                                *p(obj_gradient, x, y, w, h) += alpha * (-(2. * g_xx + g_xy + g_yx + 2. *  g_yy) / g2_norm);
                                if(x > 0) {
                                        *p(obj_gradient, x-1, y, w, h) += alpha * ((g_yx + g_xx) / g2_norm);
                                }
                                if(x < w-1) {
                                        *p(obj_gradient, x+1, y, w, h) += alpha * ((g_xx + g_xy) / g2_norm);
                                }
                                if(y > 0) {
                                        *p(obj_gradient, x, y-1, w, h) += alpha * ((g_yy + g_xy) / g2_norm);
                                }
                                if(y < h-1) {
                                        *p(obj_gradient, x, y+1, w, h) += alpha * ((g_yy + g_yx) / g2_norm);
                                }
                                if(x < w-1 && y > 0) {
                                        *p(obj_gradient, x+1, y-1, w, h) += alpha * ((-g_xy) / g2_norm);
                                }
                                if(x > 0 && y < h-1) {
                                        *p(obj_gradient, x-1, y+1, w, h) += alpha * ((-g_yx) / g2_norm);
                                }
                        }
                }
        }
        return tv2;
}

static double compute_step(
        unsigned w, unsigned h, float *in, float *out,
        float step_size, float weight, float pweight,
        int16_t *data, uint16_t quant_table[64], float *cos,
        float *temp[3],
        struct logger *log)
{
        float alpha = weight / sqrt(4. / 2.);
        float p_alpha = pweight * 2. * 255. * sqrt(2.);
        float *obj_gradient = temp[0];
        float *in_x = temp[1];
        float *in_y = temp[2];

        for(unsigned i = 0; i < h * w; i++) {
                obj_gradient[i] = 0.;
        }

        double prob_dist = pweight == 0. ? 0. : compute_step_prob(w, h, data, p_alpha, quant_table, cos, obj_gradient);

        double tv = compute_step_tv(w, h, in, obj_gradient, in_x, in_y);
#ifdef SIMD_VERIFY
        verify_compute_step_tv(w, h, tv, in, obj_gradient, in_x, in_y);
#endif

        double tv2 = alpha == 0. ? 0. : compute_step_tv2(w, h, obj_gradient, in_x, in_y, alpha);

        float norm = 0.;
        for(unsigned i = 0; i < h * w; i++) {
                norm += sqr(obj_gradient[i]);
        }
        norm = sqrt(norm);

        for(unsigned i = 0; i < h * w; i++) {
                out[i] = in[i] - step_size * (obj_gradient[i] /  norm);
        }

        double objective = (tv + alpha * tv2 + p_alpha * prob_dist) / (1. + alpha + p_alpha);
        logger_log(log, objective, prob_dist, tv, tv2);

        return objective;
}

struct compute_aux {
        float *q_min;
        float *q_max;
        float *cos;
        float *temp[3];
        float *fista;
};

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

        for(unsigned i = 0; i < 3; i++) {
                float *t = alloc_real(h * w);
                if(!t) { die("allocation error"); }
                aux->temp[i] = t;
        }

        float *fista = alloc_real(h * w);
        if(!fista) { die("allocation error"); }
        memcpy(fista, fdata, sizeof(float) * w * h);
        aux->fista = fista;
}

static void compute_aux_destroy(struct compute_aux *aux) {
        free_real(aux->cos);
        free_real(aux->q_min);
        free_real(aux->q_max);
        for(unsigned i = 0; i < 3; i++) {
                free_real(aux->temp[i]);
        }
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

void compute(struct coef *coef, struct logger *log, struct progressbar *pb, uint16_t quant_table[64], float weight, float pweight, unsigned iterations) {
        unsigned h = coef->h;
        unsigned w = coef->w;
        float *fdata = coef->fdata;
        ASSUME_ALIGNED(fdata);

        struct compute_aux aux;
        compute_aux_init(w, h, coef->data, quant_table, fdata, &aux);

        float radius = sqrt(w*h) / 2;
        for(unsigned i = 0; i < iterations; i++) {
                log->iteration = i;

                float k = i;
                for(unsigned j = 0; j < w * h; j++) {
                        aux.fista[j] = fdata[j] + (k - 2.)/(k+1.) * (fdata[j] - aux.fista[j]);
                }

                float *t = fdata;
                fdata = aux.fista;
                aux.fista = t;

                compute_step(w, h, fdata, fdata, radius / sqrt(1 + iterations), weight, pweight, coef->data, quant_table, aux.cos, aux.temp, log);
                compute_projection(w, h, fdata, aux.temp[0], aux.cos, aux.q_min, aux.q_max);

                if(pb) {
#ifdef USE_OPENMP
    #pragma omp critical(progressbar)
#endif
                        progressbar_inc(pb);
                }
        }

        coef->fdata = fdata;
        compute_aux_destroy(&aux);
}
