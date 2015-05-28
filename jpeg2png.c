#include <stdio.h>
#include <stdlib.h>
#include <stdnoreturn.h>
#include <stdbool.h>

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

static void print_quant_table(JQUANT_TBL *t) {
        for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                        printf("%d ", t->quantval[i*8+j]);
                }
                printf("\n");
        }
}

static void print_block(JCOEF b[64]) {
        for(int i = 0; i < 8; i++) {
                for(int j = 0; j < 8; j++) {
                        printf("%d ", b[i*8+j]);
                }
                printf("\n");
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

        struct jpeg_decompress_struct d;
        struct jpeg_error_mgr jerr;
        d.err = jpeg_std_error(&jerr);
        jpeg_create_decompress(&d);
        jpeg_stdio_src(&d, in);
        jpeg_read_header(&d, true);
        for(int i = 0; i < NUM_QUANT_TBLS; i++) {
                JQUANT_TBL *t = d.quant_tbl_ptrs[i];
                if(t) {
                        printf("Quant table %d\n", i);
                        print_quant_table(t);
                }
        }
        jvirt_barray_ptr *coefs = jpeg_read_coefficients(&d);
        for(int c = 0; c < d.num_components; c++) {
                jpeg_component_info *i = &d.comp_info[c];
                printf("Component %d has %d blocks y and %d blocks x\n", c, i->height_in_blocks, i->width_in_blocks);
                for(int y = 0; y < i->height_in_blocks; y+=i->v_samp_factor) {
                        for(int x = 0; x < i->width_in_blocks; x+=i->h_samp_factor) {
                                JBLOCKARRAY b = d.mem->access_virt_barray((void*)&d, coefs[c], y, i->v_samp_factor, false);
                                printf("Block %d %d %d\n", c, x, y);
                                print_block(b[0][x]);
                        }
                }
        }

        int w = d.image_width;
        int h = d.image_height;
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
        fclose(out);
        fclose(in);

        jpeg_destroy_decompress(&d);
}
