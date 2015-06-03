#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <fftw3.h>

#include "jpeg2png.h"
#include "utils.h"
#include "jpeg.h"
#include "png.h"

static void element_product(const int16_t data[64], const uint16_t quant_table[64], float *out) {
        for(int i = 0; i < 64; i++) {
                out[i] = data[i] * quant_table[i];
        }
}

static float a(int n) {
        if(n == 0) {
                return 1./sqrt(2.);
        } else {
                return 1.;
        }
}

static inline void check(int x, int y, int w, int h) {
        assert(0 <= x);
        assert(x < w);
        assert(0 <= y);
        assert(y < h);
        (void) x;
        (void) y;
        (void) w;
        (void) h;
}

static inline float *p(float *in, int x, int y, int w, int h) {
        check(x, y, w, h);
        return &in[y * w + x];
}

static void unbox(float *restrict in, float *restrict out, int w, int h) {
        assert((w & 7) == 0);
        assert((h & 7) == 0);
        for(int block_y = 0; block_y < h / 8; block_y++) {
                for(int block_x = 0; block_x < w / 8; block_x++) {
                        for(int in_y = 0; in_y < 8; in_y++) {
                                for(int in_x = 0; in_x < 8; in_x++) {
                                        *p(out, block_x * 8 + in_x, block_y * 8 + in_y, w, h) = *in++;
                                }
                        }
                }
        }
}

static void box(float *restrict in, float *restrict out, int w, int h) {
        assert((w & 7) == 0);
        assert((h & 7) == 0);
        for(int block_y = 0; block_y < h / 8; block_y++) {
                for(int block_x = 0; block_x < w / 8; block_x++) {
                        for(int in_y = 0; in_y < 8; in_y++) {
                                for(int in_x = 0; in_x < 8; in_x++) {
                                        *out++ = *p(in, block_x * 8 + in_x, block_y * 8 + in_y, w, h);
                                }
                        }
                }
        }
}

void decode_coefficients(struct coef *coef, uint16_t *quant_table) {
        coef->fdata = fftwf_alloc_real(coef->h * coef->w);
        if(!coef->fdata) { die("allocation error"); }
        int blocks = (coef->h / 8) * (coef->w / 8);
        for(int i = 0; i < blocks; i++) {
                element_product(&coef->data[i*64], quant_table, &coef->fdata[i*64]);
                for(int v = 0; v < 8; v++) {
                        for(int u = 0; u < 8; u++) {
                                coef->fdata[i*64 + v*8+u] /= a(u) * a(v);
                        }
                }
        }
        fftwf_plan p = fftwf_plan_many_r2r(
                2, (int[]){8, 8}, blocks,
                coef->fdata, (int[]){8, 8}, 1, 64,
                coef->fdata, (int[]){8, 8}, 1, 64,
                (void*)(int[]){FFTW_REDFT01, FFTW_REDFT01}, FFTW_ESTIMATE);
        fftwf_execute(p);
        fftwf_destroy_plan(p);
        for(int i = 0; i < blocks * 64; i++) {
                coef->fdata[i] /= 16.;
        }
}

static float compute_step(struct coef *coef, float weight, float step_size) {
        int w = coef->w;
        int h = coef->h;
        float *fdata = coef->fdata;
        float alpha = weight / sqrt(4. / 2.);

        float *fdata_x = NULL;
        float *fdata_y = NULL;

        if(alpha != 0.) {
                fdata_x = fftwf_alloc_real(h * w);
                fdata_y = fftwf_alloc_real(h * w);
                if(!fdata_x) { die("allocation error"); }
                if(!fdata_y) { die("allocation error"); }
        }

        float *objective_gradient = fftwf_alloc_real(h * w);
        if(!objective_gradient) { die("allocation error"); }
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

        fftwf_free(objective_gradient);
        if(alpha != 0.) {
                fftwf_free(fdata_x);
                fftwf_free(fdata_y);
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

static void compute(struct coef *coef, uint16_t quant_table[64]) {
        struct compute_projection_aux cpa;
        compute_projection_init(coef, quant_table, &cpa);

        const int N = 100;
        float tv;
        for(int i = 0; i < N; i++) {
                compute_projection(coef, &cpa);
                tv = compute_step(coef, 0.3, 1. / sqrt(1 + N));
        }

        printf("objective = %f, %f\n", tv, tv / (coef->w * coef->h) / sqrt(2.));

        compute_projection_destroy(&cpa);
}

int main(int argc, char *argv[]) {
        if(argc != 3) {
                die("usage: jpeg2png in.jpg out.png");
        }

        FILE *in = fopen(argv[1], "rb");
        if(!in) { die(NULL); }
        FILE *out = fopen(argv[2], "wb");
        if(!out) { die(NULL); }

        puts("reading jpeg");
        struct jpeg jpeg;
        read_jpeg(in, &jpeg);
        fclose(in);

        START_TIMER(decoding_coefficients);
        for(int c = 0; c < 3; c++) {
                struct coef *coef = &jpeg.coefs[c];
                decode_coefficients(coef, jpeg.quant_table[c]);
        }
        STOP_TIMER(decoding_coefficients);

        START_TIMER(unboxing);
        for(int i = 0; i < 3; i++) {
                struct coef *coef = &jpeg.coefs[i];
                float *temp = fftwf_alloc_real(coef->h * coef->w);
                if(!temp) { die("allocation error"); }

                unbox(coef->fdata, temp, coef->w, coef->h);

                fftwf_free(coef->fdata);
                coef->fdata = temp;
        }
        STOP_TIMER(unboxing);

        START_TIMER(computing);
        for(int i = 0; i < 3; i++) {
                START_TIMER(compute_1);
                struct coef *coef = &jpeg.coefs[i];
                uint16_t *quant_table = jpeg.quant_table[i];
                compute(coef, quant_table);
                STOP_TIMER(compute_1);
        }
        STOP_TIMER(computing);

        START_TIMER(adjusting_luma);
        struct coef *coef = &jpeg.coefs[0];
        for(int i = 0; i < coef->h * coef->w; i++) {
                coef->fdata[i] += 128.;
        }
        STOP_TIMER(adjusting_luma);

        START_TIMER(upsampling);
        for(int i = 0; i < 3; i++) {
                struct coef *coef = &jpeg.coefs[i];
                if(coef->h < jpeg.h) {
                        assert(coef->h / 8 == ((jpeg.h + 7) / 8 + 1) / 2);
                        float *new = fftwf_alloc_real(coef->w * coef->h * 2);
                        for(int y = 0; y < coef->h; y++) {
                                for(int x = 0; x < coef->w; x++) {
                                        new[(y*2)*coef->w + x] = coef->fdata[y * coef->w + x];
                                }
                                for(int x = 0; x < coef->w; x++) {
                                        if(y == coef->h - 1) {
                                                new[(y*2+1)*coef->w + x] = coef->fdata[y * coef->w + x];
                                        } else {
                                                new[(y*2+1)*coef->w + x] = (coef->fdata[y * coef->w + x] + coef->fdata[(y+1) * coef->w + x]) / 2;
                                        }
                                }
                        }
                        fftwf_free(coef->fdata);
                        coef->fdata = new;
                        coef->h = coef->h * 2;
                }
                if(coef->w < jpeg.w) {
                        assert(coef->w / 8 == ((jpeg.w + 7) / 8 + 1) / 2);
                        float *new = fftwf_alloc_real(coef->w * coef->h * 2);
                        for(int y = 0; y < coef->h; y++) {
                                for(int x = 0; x < coef->w; x++) {
                                        new[(y*coef->w + x)*2] = coef->fdata[y * coef->w + x];
                                        if(x == coef->w - 1) {
                                                new[(y*coef->w + x)*2 + 1] = coef->fdata[y * coef->w + x];
                                        } else {
                                                new[(y*coef->w + x)*2 + 1] = (coef->fdata[y * coef->w + x] + coef->fdata[y * coef->w + x + 1]) / 2;
                                        }
                                }
                        }
                        fftwf_free(coef->fdata);
                        coef->fdata = new;
                        coef->w = coef->w * 2;
                }
        }
        STOP_TIMER(upsampling);

        START_TIMER(writing_png);
        write_png(out, jpeg.w, jpeg.h, &jpeg.coefs[0], &jpeg.coefs[1], &jpeg.coefs[2]);
        fclose(out);
        STOP_TIMER(writing_png);

        for(int i = 0; i < 3; i++) {
                fftwf_free(jpeg.coefs[i].fdata);
                free(jpeg.coefs[i].data);
        }
}
