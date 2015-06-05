#include <stdint.h>
#include <string.h>
#include <fftw3.h>

#include "jpeg2png.h"
#include "compute.h"
#include "utils.h"
#include "box.h"
#include "logger.h"


static float compute_step_tv(unsigned w, unsigned h, float *in, float *objective_gradient, float *in_x, float *in_y) {
        float tv = 0.;
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
                                *p(objective_gradient, x, y, w, h) += -(g_x + g_y) / g_norm;
                                if(x < w-1) {
                                        *p(objective_gradient, x+1, y, w, h) += g_x / g_norm;
                                }
                                if(y < h-1) {
                                        *p(objective_gradient, x, y+1, w, h) += g_y / g_norm;
                                }
                        }
                        *p(in_x, x, y, w, h) = g_x;
                        *p(in_y, x, y, w, h) = g_y;
                }
        }
        return tv;
}

static float compute_step_tv2(unsigned w, unsigned h, float *objective_gradient, float *in_x, float *in_y, float alpha) {
        float tv2 = 0.;
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
                                *p(objective_gradient, x, y, w, h) += alpha * (-(2. * g_xx + g_xy + g_yx + 2. *  g_yy) / g2_norm);
                                if(x > 0) {
                                        *p(objective_gradient, x-1, y, w, h) += alpha * ((g_yx + g_xx) / g2_norm);
                                }
                                if(x < w-1) {
                                        *p(objective_gradient, x+1, y, w, h) += alpha * ((g_xx + g_xy) / g2_norm);
                                }
                                if(y > 0) {
                                        *p(objective_gradient, x, y-1, w, h) += alpha * ((g_yy + g_xy) / g2_norm);
                                }
                                if(y < h-1) {
                                        *p(objective_gradient, x, y+1, w, h) += alpha * ((g_yy + g_yx) / g2_norm);
                                }
                                if(x < w-1 && y > 0) {
                                        *p(objective_gradient, x+1, y-1, w, h) += alpha * ((-g_xy) / g2_norm);
                                }
                                if(x > 0 && y < h-1) {
                                        *p(objective_gradient, x-1, y+1, w, h) += alpha * ((-g_yx) / g2_norm);
                                }
                        }
                }
        }
        return tv2;
}

static float compute_step(unsigned w, unsigned h, float *in, float *out, float step_size, float weight, float *objective_gradient, float *in_x, float *in_y, struct logger *log) {
        float alpha = weight / sqrt(4. / 2.);
        for(unsigned i = 0; i < h * w; i++) {
                objective_gradient[i] = 0.;
        }

        float tv = compute_step_tv(w, h, in, objective_gradient, in_x, in_y);

        float tv2 = alpha == 0. ? 0. : compute_step_tv2(w, h, objective_gradient, in_x, in_y, alpha);

        float norm = 0.;
        for(unsigned i = 0; i < h * w; i++) {
                norm += sqr(objective_gradient[i]);
        }
        norm = sqrt(norm);

        for(unsigned i = 0; i < h * w; i++) {
                out[i] = in[i] - step_size * (objective_gradient[i] /  norm);
        }

        float objective = (tv + alpha * tv2) / (alpha + 1.);
        logger_log(log, objective, tv, tv2);

        return objective;
}

struct compute_projection_aux {
        float *q_min;
        float *q_max;
        float *temp;
        fftwf_plan dct;
        fftwf_plan idct;
};

static void compute_projection_init(unsigned w, unsigned h, int16_t *data, uint16_t quant_table[64], struct compute_projection_aux *aux) {
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

        for(unsigned i = 0; i < blocks; i++) {
                for(unsigned v = 0; v < 8; v++) {
                        for(unsigned u = 0; u < 8; u++) {
                                q_max[i*64 + v*8+u] /= a(u) * a(v);
                                q_min[i*64 + v*8+u] /= a(u) * a(v);
                        }
                }
        }

        aux->q_min = q_min;
        aux->q_max = q_max;

        float *temp = alloc_real(h * w);
        if(!temp) { die("allocation error"); }

        aux->temp = temp;

        fftwf_plan dct = fftwf_plan_many_r2r(
                2, (int[]){8, 8}, blocks,
                temp, (int[]){8, 8}, 1, 64,
                temp, (int[]){8, 8}, 1, 64,
                (void*)(int[]){FFTW_REDFT10, FFTW_REDFT10}, FFTW_ESTIMATE);

        aux->dct = dct;

        fftwf_plan idct = fftwf_plan_many_r2r(
                2, (int[]){8, 8}, blocks,
                temp, (int[]){8, 8}, 1, 64,
                temp, (int[]){8, 8}, 1, 64,
                (void*)(int[]){FFTW_REDFT01, FFTW_REDFT01}, FFTW_ESTIMATE);

        aux->idct = idct;
}

static void compute_projection_destroy(struct compute_projection_aux *aux) {
        fftwf_destroy_plan(aux->idct);
        fftwf_destroy_plan(aux->dct);
        free(aux->temp);
        free(aux->q_min);
        free(aux->q_max);
}

static void compute_projection(unsigned w, unsigned h, float *fdata, struct compute_projection_aux *aux) {
        float *temp = aux->temp;

        unsigned blocks = (h / 8) * (w / 8);

        box(fdata, temp, w, h);

        fftwf_execute(aux->dct);
        for(unsigned i = 0; i < h * w; i++) {
                temp[i] /= 16.;
        }

        for(unsigned i = 0; i < h * w; i++) {
                temp[i] = fmin(fmax(temp[i], aux->q_min[i]), aux->q_max[i]);
        }

        fftwf_execute(aux->idct);
        for(unsigned i = 0; i < blocks * 64; i++) {
                temp[i] /= 16.;
        }

        unbox(temp, fdata, w, h);
}

void compute(struct coef *coef, struct logger *log, uint16_t quant_table[64], float weight, unsigned iterations) {
        unsigned h = coef->h;
        unsigned w = coef->w;
        float *fdata = coef->fdata;
        ASSUME_ALIGNED(fdata);

        struct compute_projection_aux cpa;
        compute_projection_init(w, h, coef->data, quant_table, &cpa);

        float *temp_x = alloc_real(h * w);
        if(!temp_x) { die("allocation error"); }
        float *temp_y = alloc_real(h * w);
        if(!temp_y) { die("allocation error"); }
        float *temp_gradient = alloc_real(h * w);
        if(!temp_gradient) { die("allocation error"); }

        float *temp_fista = alloc_real(h * w);
        if(!temp_fista) { die("allocation error"); }
        memcpy(temp_fista, fdata, sizeof(float) * w * h);

        float radius = sqrt(w*h) / 2;
        for(unsigned i = 0; i < iterations; i++) {
                log->iteration = i;

                float k = i;
                for(unsigned j = 0; j < w * h; j++) {
                        temp_fista[j] = fdata[j] + (k - 2.)/(k+1.) * (fdata[j] - temp_fista[j]);
                }

                compute_step(w, h, temp_fista, temp_fista, radius / sqrt(1 + iterations), weight, temp_x, temp_y, temp_gradient, log);
                compute_projection(w, h, temp_fista, &cpa);

                float *t = fdata;
                fdata = temp_fista;
                temp_fista = t;
        }

        coef->fdata = fdata;
        free(temp_x);
        free(temp_y);
        free(temp_gradient);
        free(temp_fista);

        compute_projection_destroy(&cpa);
}
