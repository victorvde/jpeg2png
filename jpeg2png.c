#include <stdio.h>
#include <stdlib.h>
#include <stdnoreturn.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include <jpeglib.h>
#include <png.h>

static noreturn void die(const char *msg)  {
        if(msg) {
                fprintf(stderr, "jpeg2png: %s\n");
        } else {
                perror("jpeg2png");
        }
        exit(EXIT_FAILURE);
}

static void print_quant_table(uint16_t *t) {
        for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                        printf("%d ", t[i*8+j]);
                }
                printf("\n");
        }
}

static void print_block(int16_t b[64]) {
        for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                        printf("%d ", b[i*8+j]);
                }
                printf("\n");
        }
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
                memcpy(&jpeg->quant_table[c], t->quantval, sizeof(uint16_t) * 64);
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
                for(int y = 0; y < i->height_in_blocks; y+=i->v_samp_factor) {
                        for(int x = 0; x < i->width_in_blocks; x+=i->h_samp_factor) {
                                JBLOCKARRAY b = d.mem->access_virt_barray((void*)&d, coefs[c], y, i->v_samp_factor, false);
                                memcpy(data, b[0][x], 64 * sizeof(int16_t));
                                data += 64;
                        }
                }
        }
        jpeg_destroy_decompress(&d);
}

static void write_png(FILE *out, int w, int h) {
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
        png_byte *rows[h];
        for(int i = 0; i < h; i++) {
                rows[i] = &image_data[w * 3];
        }
        png_write_image(png_ptr, rows);
        png_write_end(png_ptr, info_ptr);
        free(image_data);
        png_destroy_write_struct(&png_ptr, &info_ptr);
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
        for(int c = 0; c < 3; c++) {
                printf("Quant table %d\n", c);
                print_quant_table(jpeg.quant_table[c]);
        }

        for(int c = 0; c < 3; c++) {
                struct coefs coefs = jpeg.coefs[c];
                int16_t *d = coefs.data;
                for(int y = 0; y < coefs.h; y += 8) {
                        for(int x = 0; x < coefs.w; x += 8) {
                                printf("Block at %dx%d\n", x, y);
                                print_block(d);
                                d += 64;
                        }
                }
        }

        write_png(out, jpeg.w, jpeg.h);
        fclose(out);
}
