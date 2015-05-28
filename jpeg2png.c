#include <stdio.h>
#include <stdlib.h>
#include <stdnoreturn.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>

#include <jpeglib.h>
#include <png.h>

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

static int clamp(int x) {
        if(x < 0) {
                return 0;
        } else if(x > 255) {
                return 255;
        } else {
                return x;
        }
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

        png_byte *rows[h];
        for(int i = 0; i < h; i++) {
                rows[i] = &image_data[i * w * 3];
        }
        png_write_image(png_ptr, rows);
        png_write_end(png_ptr, info_ptr);
        free(image_data);
        png_destroy_write_struct(&png_ptr, &info_ptr);
}

static void element_product(int16_t data[64], const uint16_t quant_table[64]) {
        for(int i = 0; i < 64; i++) {
                data[i] *= quant_table[i];
        }
}

static float a(int n) {
        if(n == 0) {
                return 1./sqrt(2.);
        } else {
                return 1.;
        }
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
static void idct(int16_t data[64], float out[64]) {
        // super slow idct
        for(int y = 0; y < 8; y++) {
                for(int x = 0; x < 8; x++) {
                        float sum = 0.;
                        for(int u = 0; u < 8; u++) {
                                for(int v = 0; v < 8; v++) {
                                        sum += a(u) * a(v) * data[v*8+u] * cos((2*x+1)*u*M_PI / 16.) * cos((2*y+1)*v*M_PI / 16.);
                                }
                        }
                        out[y*8+x] = sum / 4.;
                }
        }
}

static void unbox(const float *in, float *out, int w, int h) {
        for(int y = 0; y < h; y+=8) {
                for(int x = 0; x < w; x+=8) {
                        for(int iy = 0; iy < 8; iy++) {
                                for(int ix = 0; ix < 8; ix++) {
                                        out[(y+iy)*w+x+ix] = *in++;
                                }
                        }
                }
        }
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

        int w = jpeg.w;
        int h = jpeg.h;

        puts("decoding coefficients");
        for(int c = 0; c < 3; c++) {
                struct coef *coef = &jpeg.coefs[c];
                if(coef->h < h) { die("channel %d too short! %d", c, coef->h); }
                if(coef->w < w) { die("channel %d too thin! %d", c, coef->w); }
                coef->fdata = calloc(coef->h * coef->w, sizeof(float));
                if(!coef->fdata) { die("allocation error"); }
                int16_t *d = coef->data;
                float *f = coef->fdata;
                for(int y = 0; y < coef->h; y += 8) {
                        for(int x = 0; x < coef->w; x += 8) {
                                element_product(d, jpeg.quant_table[c]);
                                idct(d, f);
                                d += 64;
                                f += 64;
                        }
                }
        }

        puts("unboxing");
        for(int i = 0; i < 3; i++) {
                struct coef *coef = &jpeg.coefs[i];
                float *temp = calloc(sizeof(float), coef->h * coef->w);
                if(!temp) { die("allocation error"); }

                unbox(coef->fdata, temp, w, h);
                free(coef->fdata);
                coef->fdata = temp;
        }

        puts("adjusting luma");
        for(int i = 0; i < h * w; i++) {
                jpeg.coefs[0].fdata[i] += 128.;
        }

        puts("writing png");
        struct coef yc = jpeg.coefs[0];
        struct coef cbc = jpeg.coefs[1];
        struct coef crc = jpeg.coefs[2];
        write_png(out, w, h, yc.fdata, cbc.fdata, crc.fdata);
        fclose(out);

        for(int i = 0; i < 3; i++) {
                free(jpeg.coefs[i].fdata);
        }
}
