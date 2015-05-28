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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static noreturn void die(const char *msg, ...)  {
        if(msg) {
                fprintf(stderr, "jpeg2png: ");
                va_list l;
                va_start(l, msg);
                vfprintf(stderr, msg, l);
                va_end(l);
        } else {
                perror("jpeg2png");
        }
        exit(EXIT_FAILURE);
}

struct coefs {
        int h;
        int w;
        int16_t *data;
};

struct jpeg {
        int h;
        int w;
        uint16_t quant_table[3][64];
        struct coefs coefs[3];
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

static void write_png(FILE *out, int w, int h, float *fdata[3]) {
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
        if(!image_data) { die("Could not allocate image data");}

        for(int y = 0; y < h; y++) {
                for(int x = 0; x < w; x++) {
                        int i = y * w + x;
                        image_data[i*3] = clamp(fdata[0][i] + 1.402 * fdata[2][i]);
                        image_data[i*3+1] = clamp(fdata[0][i] - 0.34414 * fdata[1][i] - 0.71414 * fdata[2][i]);
                        image_data[i*3+2] = clamp(fdata[0][i] + 1.772 * fdata[1][i]);
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

static void idct(int16_t data[64], int w, int h, int bx, int by, float *out) {
        // super slow
        for(int x = 0; x < 8; x++) {
                for(int y = 0; y < 8; y++) {
                        int nx = bx + x;
                        int ny = by + y;
                        if(nx < w && ny < h) {
                                float sum = 0.;
                                for(int u = 0; u < 8; u++) {
                                        for(int v = 0; v < 8; v++) {
                                                sum += a(u) * a(v) * data[v*8+u] * cos((2*x+1)*u*M_PI / 16.) * cos((2*y+1)*v*M_PI / 16.);
                                        }
                                }
                                out[ny*w+nx] = sum / 4.;
                        }
                }
        }
}

int main(int argc, char *argv[]) {
        if(argc != 3) {
                fprintf(stderr, "usage: jpeg2png in.jpg out.png\n");
                exit(EXIT_FAILURE);
        }

        FILE *in = fopen(argv[1], "rb");
        if(!in) { die(NULL); }
        FILE *out = fopen(argv[2], "wb");
        if(!out) { die(NULL); }

        struct jpeg jpeg;
        read_jpeg(in, &jpeg);
        fclose(in);

        float *fdata[3];
        for(int i = 0; i < 3; i++) {
                fdata[i] = calloc(sizeof(float), jpeg.h * jpeg.w);
        }
        for(int c = 0; c < 3; c++) {
                struct coefs coefs = jpeg.coefs[c];
                int16_t *d = coefs.data;
                if(coefs.h < jpeg.h) { die("channel %d too short! %d", c, coefs.h); }
                if(coefs.w < jpeg.w) { die("channel %d too thin! %d", c, coefs.w); }
                for(int y = 0; y < coefs.h; y += 8) {
                        for(int x = 0; x < coefs.w; x += 8) {
                                element_product(d, jpeg.quant_table[c]);
                                idct(d, jpeg.w, jpeg.h, x, y, fdata[c]);
                                d += 64;
                        }
                }
        }
        for(int i = 0; i < jpeg.h * jpeg.w; i++) {
                fdata[0][i] += 128.;
        }

        write_png(out, jpeg.w, jpeg.h, fdata);
        fclose(out);
}
