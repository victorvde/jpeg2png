#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include <jpeglib.h>

#include "jpeg.h"
#include "utils.h"
#include "ooura/dct.h"

// clean up progress bar when printing warnings
static void die_output_message(struct jpeg_common_struct *c) {
        die_message_start();
        char error_message[JMSG_LENGTH_MAX];
        c->err->format_message(c, error_message);
        fprintf(stderr, "libjpeg error: %s\n", error_message);
}

// read JPEG file DCT coefficients and quantization tables
void read_jpeg(FILE *in, struct jpeg *jpeg) {
        struct jpeg_decompress_struct d;
        struct jpeg_error_mgr jerr;
        d.err = jpeg_std_error(&jerr);
        d.err->output_message = die_output_message;
        jpeg_create_decompress(&d);
        jpeg_stdio_src(&d, in);
        jpeg_read_header(&d, true);

        jpeg->h = d.image_height;
        jpeg->w = d.image_width;

        if(d.num_components != 3) { die("only 3 component jpegs are supported"); }

        for(int c = 0; c < d.num_components; c++) {
                unsigned i = d.comp_info[c].quant_tbl_no;
                if(i >= NUM_QUANT_TBLS) { die("weird jpeg: invalid quant_tbl_no"); }
                JQUANT_TBL *t = d.quant_tbl_ptrs[i];
                if(!t) { die("weird jpeg: no quant table pointer"); }
                for(unsigned j = 0; j < 64; j++) {
                        if(t->quantval[j] == 0) {
                                die("invalid quantization table");
                        }
                }
                memcpy(&(jpeg->coefs[c].quant_table), t->quantval, sizeof(uint16_t) * 64);
        }

        jvirt_barray_ptr *coefs = jpeg_read_coefficients(&d);
        for(int c = 0; c < d.num_components; c++) {
                jpeg_component_info *i = &d.comp_info[c];
                unsigned h = i->height_in_blocks * 8;
                unsigned w = i->width_in_blocks * 8;
                struct coef *coef = &jpeg->coefs[c];
                coef->w = w;
                coef->h = h;
                coef->w_samp = d.max_h_samp_factor / i->h_samp_factor;
                coef->h_samp = d.max_v_samp_factor / i->v_samp_factor;
                if(coef->h / 8 != (jpeg->h / coef->h_samp + 7) / 8) {
                        die("jpeg invalid coef h size");
                }
                if(coef->w / 8 != (jpeg->w / coef->w_samp + 7) / 8) {
                        die("jpeg invalid coef w size");
                }
                if(SIZE_MAX / coef->h / coef->w / coef->h_samp / coef->w_samp < 6) {
                        die("jpeg is too big to fit in memory");
                }
                int16_t *data = malloc(sizeof(*data) * h * w);
                if(!data) { die("could not allocate memory for coefs"); }
                coef->data = data;
                for(unsigned y = 0; y < h / 8; y++) {
                        JBLOCKARRAY b = d.mem->access_virt_barray((void*)&d, coefs[c], y, 1, false);
                        for(unsigned x = 0; x < w / 8; x++) {
                                memcpy(data, b[0][x], 64 * sizeof(int16_t));
                                data += 64;
                        }
                }
        }
        jpeg_destroy_decompress(&d);
}

// decode DCT coefficients into image data
void decode_coefficients(struct coef *coef) {
        coef->fdata = alloc_simd(sizeof(float) * coef->h * coef->w);
        unsigned blocks = (coef->h / 8) * (coef->w / 8);
        for(unsigned i = 0; i < blocks; i++) {
                for(unsigned j = 0; j < 64; j++) {
                        coef->fdata[i*64+j] = coef->data[i*64+j] * coef->quant_table[j];
                }
                idct8x8s(&(coef->fdata[i*64]));
        }
}
