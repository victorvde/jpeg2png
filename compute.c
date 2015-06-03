#include <stdint.h>
#include <fftw3.h>

#include "jpeg2png.h"
#include "compute.h"
#include "utils.h"
#include "box.h"

static float compute_step(struct coef *coef, float *objective_gradient, float *fdata_x, float *fdata_y, float weight, float step_size) {
        int w = coef->w;
        int h = coef->h;
        float *fdata = coef->fdata;
        float alpha = weight / sqrt(4. / 2.);

        for(int i = 0; i < h * w; i++) {
                objective_gradient[i] = 0.;
        }

        float tv = 0.;
        for(int y = 0; y < h; y++) {
                for(int x = 0; x < w; x++) {
                        // forward gradient x
                        float g_x = x >= w-1 ? 0. : *p(fdata, x+1, y, w, h) - *p(fdata, x, y, w, h);
                        // forward gradient y
                        float g_y = y >= h-1 ? 0. : *p(fdata, x, y+1, w, h) - *p(fdata, x, y, w, h);
                        // norm
                        float g_norm = sqrt(g_x * g_x + g_y * g_y);
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
                        if(alpha != 0.) {
                                *p(fdata_x, x, y, w, h) = g_x;
                                *p(fdata_y, x, y, w, h) = g_y;
                        }
                }
        }
        // printf("tv = %f, %f\n", tv, tv / (w * h) / sqrt(2.));

        float tv2 = 0;
        if(alpha != 0.) {
                for(int y = 0; y < h; y++) {
                        for(int x = 0; x < w; x++) {
                                // backward x
                                float g_xx = x <= 0 ? 0. : *p(fdata_x, x, y, w, h) - *p(fdata_x, x-1, y, w, h);
                                // backward x
                                float g_yx = x <= 0 ? 0. : *p(fdata_y, x, y, w, h) - *p(fdata_y, x-1, y, w, h);
                                // backward y
                                float g_xy = y <= 0 ? 0. : *p(fdata_x, x, y, w, h) - *p(fdata_x, x, y-1, w, h);
                                // backward y
                                float g_yy = y <= 0 ? 0. : *p(fdata_y, x, y, w, h) - *p(fdata_y, x, y-1, w, h);
                                // norm
                                float g2_norm = sqrt(g_xx * g_xx + g_yx * g_yx + g_xy * g_xy + g_yy * g_yy);
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
                // printf("tv2 = %f, %f\n", tv2, tv2 / (w * h));
        }

        for(int i = 0; i < h * w; i++) {
                fdata[i] -= step_size * (objective_gradient[i] / (alpha + 1.));
        }

        return (tv + alpha * tv2) / (alpha + 1.);
}

struct compute_projection_aux {
        float * restrict q_min;
        float * restrict q_max;
        float *temp;
        fftwf_plan dct;
        fftwf_plan idct;
};

static void compute_projection_init(struct coef *coef, uint16_t quant_table[64],struct compute_projection_aux *aux) {
        int w = coef->w;
        int h = coef->h;

        float *q_max = fftwf_alloc_real(h * w);
        if(!q_max) { die("allocation error"); }
        float *q_min = fftwf_alloc_real(h * w);
        if(!q_min) { die("allocation error"); }
        int blocks = (h / 8) * (w / 8);

        for(int i = 0; i < blocks; i++) {
                for(int j = 0; j < 64; j++) {
                       q_max[i*64+j] = (coef->data[i*64+j] + 0.5) * quant_table[j];
                       q_min[i*64+j] = (coef->data[i*64+j] - 0.5) * quant_table[j];
                }
        }

        for(int i = 0; i < blocks; i++) {
                for(int v = 0; v < 8; v++) {
                        for(int u = 0; u < 8; u++) {
                                q_max[i*64 + v*8+u] /= a(u) * a(v);
                                q_min[i*64 + v*8+u] /= a(u) * a(v);
                        }
                }
        }

        aux->q_min = q_min;
        aux->q_max = q_max;

        float *temp = fftwf_alloc_real(h * w);
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
        fftwf_free(aux->temp);
        fftwf_free(aux->q_min);
        fftwf_free(aux->q_max);
}

static void compute_projection(struct coef *coef, struct compute_projection_aux *aux) {
        int w = coef->w;
        int h = coef->h;
        float *fdata = coef->fdata;

        float *temp = aux->temp;

        int blocks = (h / 8) * (w / 8);

        box(fdata, temp, w, h);

        fftwf_execute(aux->dct);
        for(int i = 0; i < h * w; i++) {
                temp[i] /= 16.;
        }

        for(int i = 0; i < h * w; i++) {
                temp[i] = CLAMP(temp[i], aux->q_min[i], aux->q_max[i]);
        }

        fftwf_execute(aux->idct);
        for(int i = 0; i < blocks * 64; i++) {
                temp[i] /= 16.;
        }

        unbox(temp, fdata, w, h);
}

void compute(struct coef *coef, uint16_t quant_table[64], float weight, int iterations) {
        struct compute_projection_aux cpa;
        compute_projection_init(coef, quant_table, &cpa);

        int h = coef->h;
        int w = coef->w;

        float *fdata_x = fftwf_alloc_real(h * w);
        if(!fdata_x) { die("allocation error"); }
        float *fdata_y = fftwf_alloc_real(h * w);
        if(!fdata_y) { die("allocation error"); }
        float *objective_gradient = fftwf_alloc_real(h * w);
        if(!objective_gradient) { die("allocation error"); }

        float objective = 0.;
        for(int i = 0; i < iterations; i++) {
                compute_projection(coef, &cpa);
                objective = compute_step(coef, fdata_x, fdata_y, objective_gradient, weight, 1. / sqrt(1 + iterations));
        }

        printf("objective = %f, %f\n", objective, objective / (coef->w * coef->h) / sqrt(2.));

        fftwf_free(fdata_x);
        fftwf_free(fdata_y);
        fftwf_free(objective_gradient);
        compute_projection_destroy(&cpa);
}
