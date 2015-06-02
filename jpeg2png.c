#include <stdio.h>
#include <stdlib.h>
#include <stdnoreturn.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <assert.h>
#include <time.h>

#include <jpeglib.h>
#include <png.h>

#include <fftw3.h>

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))

static noreturn void die(const char *msg, ...)  {
        if(msg) {
                fprintf(stderr, "jpeg2png: ");
                va_list l;
                va_start(l, msg);
                vfprintf(stderr, msg, l);
                va_end(l);
                fprintf(stderr, "\n");
        } else {
                perror("jpeg2png");
        }
        exit(EXIT_FAILURE);
}

struct coef {
        int h;
        int w;
        int16_t *data;
        float *fdata;
};

struct jpeg {
        int h;
        int w;
        uint16_t quant_table[3][64];
        struct coef coefs[3];
};

static void read_jpeg(FILE *in, struct jpeg *jpeg) {
        struct jpeg_decompress_struct d;
        struct jpeg_error_mgr jerr;
        d.err = jpeg_std_error(&jerr);
        jpeg_create_decompress(&d);
        jpeg_stdio_src(&d, in);
        jpeg_read_header(&d, true);

        jpeg->h = d.image_height;
        jpeg->w = d.image_width;

        if(d.num_components != 3) { die("only 3 component jpegs are supported"); }

        for(int c = 0; c < d.num_components; c++) {
                int i = d.comp_info[c].quant_tbl_no;
                JQUANT_TBL *t = d.quant_tbl_ptrs[i];
                memcpy(jpeg->quant_table[c], t->quantval, sizeof(uint16_t) * 64);
        }

        jvirt_barray_ptr *coefs = jpeg_read_coefficients(&d);
        for(int c = 0; c < d.num_components; c++) {
                jpeg_component_info *i = &d.comp_info[c];
                int h = ((i->height_in_blocks + i->v_samp_factor - 1) / i->v_samp_factor) * 8;
                int w = ((i->width_in_blocks + i->h_samp_factor - 1) / i->h_samp_factor) * 8;
                int16_t *data = malloc(w * h * sizeof(int16_t));
                jpeg->coefs[c].w = w;
                jpeg->coefs[c].h = h;
                jpeg->coefs[c].data = data;
                if(!data) { die("could not allocate memory for coefs"); }
                for(unsigned y = 0; y < i->height_in_blocks; y+=i->v_samp_factor) {
                        for(unsigned x = 0; x < i->width_in_blocks; x+=i->h_samp_factor) {
                                JBLOCKARRAY b = d.mem->access_virt_barray((void*)&d, coefs[c], y, i->v_samp_factor, false);
                                memcpy(data, b[0][x], 64 * sizeof(int16_t));
                                data += 64;
                        }
                }
        }
        jpeg_destroy_decompress(&d);
}

static float clamp(float x) {
        return CLAMP(x, 0., 255.);
}

static void write_png(FILE *out, int w, int h, float *y, float *cb, float *cr) {
        png_struct *png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if(!png_ptr) { die("could not initialize PNG write struct"); }
        png_info *info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) { die("could not initialize PNG info struct"); }
        if (setjmp(png_jmpbuf(png_ptr)))
        {
                /* If we get here, we had a problem writing the file */
                die("problem while writing PNG file");
        }
        png_init_io(png_ptr, out);
        png_set_IHDR(png_ptr, info_ptr, w, h, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_write_info(png_ptr, info_ptr);
        png_byte *image_data = calloc(sizeof(png_byte), h * w * 3);
        if(!image_data) { die("could not allocate image data");}

        for(int yi = 0; yi < h; yi++) {
                for(int xi = 0; xi < w; xi++) {
                        int i = yi * w + xi;
                        image_data[i*3] = clamp(y[i] + 1.402 * cr[i]);
                        image_data[i*3+1] = clamp(y[i] - 0.34414 * cb[i] - 0.71414 * cr[i]);
                        image_data[i*3+2] = clamp(y[i] + 1.772 * cb[i]);
                }
        }

        png_byte **rows = malloc(sizeof(*rows) * h);
        if(!rows) { die("allocation failure"); }
        for(int i = 0; i < h; i++) {
                rows[i] = &image_data[i * w * 3];
        }
        png_write_image(png_ptr, rows);
        free(rows);
        png_write_end(png_ptr, info_ptr);
        free(image_data);
        png_destroy_write_struct(&png_ptr, &info_ptr);
}

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
        for(int y = 0; y < h; y++) {
                for(int x = 0; x < w; x++) {
                        int y_block = y / 8;
                        int y_in = y & 7;
                        int block_width = w / 8;
                        int x_block = x / 8;
                        int x_in = x & 7;
                        float t = in[(y_block*block_width + x_block) * 64 + y_in * 8 + x_in];
                        *p(out, x, y, w, h) = t;
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

static inline float gradient_x_forward(float *in, int x, int y, int w, int h) {
        if(x == w-1) {
                return 0.;
        } else {
                return *p(in, x+1, y, w, h) - *p(in, x, y, w, h);
        }
}
static inline float gradient_x_backward(float *in, int x, int y, int w, int h) {
        if(x == 0) {
                return 0.;
        } else {
                return *p(in, x, y, w, h) - *p(in, x-1, y, w, h);
        }
}

static inline float gradient_y_forward(float *in, int x, int y, int w, int h) {
        if(y == h-1) {
                return 0.;
        } else {
                return *p(in, x, y+1, w, h) - *p(in, x, y, w, h);
        }
}
static inline float gradient_y_backward(float *in, int x, int y, int w, int h) {
        check(x, y, w, h);
        if(y == 0) {
                return 0.;
        } else {
                return *p(in, x, y, w, h) - *p(in, x, y-1, w, h);
        }
}

#define START_TIMER(n) clock_t macro_timer_##n = start_timer(#n);
static clock_t start_timer(const char *name) {
        (void) name;
        return clock();
}

#define STOP_TIMER(n) stop_timer(macro_timer_##n, #n);
static void stop_timer(clock_t t, const char *n) {
        clock_t diff = clock() - t;
        int msec = diff * 1000 / CLOCKS_PER_SEC;
        printf("%s: %d ms\n", n, msec);
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

static float compute_step(struct coef *coef, float step_size) {
        int w = coef->w;
        int h = coef->h;
        float *fdata = coef->fdata;

        // float *fdata_x = fftwf_alloc_real(h * w);
        // if(!fdata_x) { die("allocation error"); }
        // float *fdata_y = fftwf_alloc_real(h * w);
        // if(!fdata_y) { die("allocation error"); }

        float *objective_gradient = fftwf_alloc_real(h * w);
        if(!objective_gradient) { die("allocation error"); }
        for(int i = 0; i < h * w; i++) {
                objective_gradient[i] = 0.;
        }

        float tv = 0;
        for(int y = 0; y < h; y++) {
                for(int x = 0; x < w; x++) {
                        float g_x = gradient_x_forward(fdata, x, y, w, h);
                        float g_y = gradient_y_forward(fdata, x, y, w, h);
                        float g_norm = hypotf(g_x, g_y);
                        tv += g_norm;
                        if(g_norm != 0) {
                                if(g_x != 0.) {
                                        *p(objective_gradient, x, y, w, h) += g_x / g_norm;
                                        *p(objective_gradient, x+1, y, w, h) -= g_x / g_norm;
                                }
                                if(g_y != 0.) {
                                        *p(objective_gradient, x, y, w, h) += g_y / g_norm;
                                        *p(objective_gradient, x, y+1, w, h) -= g_y / g_norm;
                                }
                        }
                        // *p(fdata_x, x, y, w, h) = g_x;
                        // *p(fdata_y, x, y, w, h) = g_y;
                }
        }
        // printf("tv = %f, %f\n", tv, tv / (w * h) / sqrt(2.));
        for(int i = 0; i < h * w; i++) {
                fdata[i] += step_size * objective_gradient[i];
        }

        // float tv2 = 0;
        // for(int y = 0; y < h; y++) {
        //         for(int x = 0; x < w; x++) {
        //                 float g_xx = gradient_x_backward(fdata_x, x, y, w, h);
        //                 float g_xy = gradient_y_backward(fdata_x, x, y, w, h);
        //                 float g_yy = gradient_y_backward(fdata_y, x, y, w, h);

        //                 tv2 += sqrt(g_xx * g_xx + 2 * g_xy * g_xy + g_yy * g_yy);
        //         }
        // }
        // printf("tv2 = %f, %f\n", tv2, tv2 / (w * h) / sqrt(4.));

        fftwf_free(objective_gradient);
        // fftwf_free(fdata_x);
        // fftwf_free(fdata_y);

        return tv;
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
                tv = compute_step(coef, 1. / sqrt(1 + N));
        }

        printf("TV = %f, %f\n", tv, tv / (coef->w * coef->h) / sqrt(2.));

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
                if(coef->h < jpeg.h) { die("channel %d too short! %d", c, coef->h); }
                if(coef->w < jpeg.w) { die("channel %d too thin! %d", c, coef->w); }
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

        START_TIMER(writing_png);
        struct coef yc = jpeg.coefs[0];
        printf("luma: %d x %d\n", yc.w, yc.h);
        printf("jpeg: %d x %d\n", jpeg.w, jpeg.h);
        struct coef cbc = jpeg.coefs[1];
        struct coef crc = jpeg.coefs[2];
        write_png(out, jpeg.w, jpeg.h, yc.fdata, cbc.fdata, crc.fdata);
        fclose(out);
        STOP_TIMER(writing_png);

        for(int i = 0; i < 3; i++) {
                fftwf_free(jpeg.coefs[i].fdata);
                free(jpeg.coefs[i].data);
        }
}
