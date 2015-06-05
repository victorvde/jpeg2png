#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

#include <jpeglib.h>
#include <fftw3.h>

#include "jpeg.h"
#include "utils.h"

void read_jpeg(FILE *in, struct jpeg *jpeg) {
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
                unsigned i = d.comp_info[c].quant_tbl_no;
                JQUANT_TBL *t = d.quant_tbl_ptrs[i];
                memcpy(jpeg->quant_table[c], t->quantval, sizeof(uint16_t) * 64);
        }

        jvirt_barray_ptr *coefs = jpeg_read_coefficients(&d);
        for(int c = 0; c < d.num_components; c++) {
                jpeg_component_info *i = &d.comp_info[c];
                unsigned h = i->height_in_blocks * 8;
                unsigned w = i->width_in_blocks * 8;
                int16_t *data = malloc(w * h * sizeof(int16_t));
                jpeg->coefs[c].w = w;
                jpeg->coefs[c].h = h;
                jpeg->coefs[c].data = data;
                if(!data) { die("could not allocate memory for coefs"); }
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

void decode_coefficients(struct coef *coef, uint16_t *quant_table) {
        coef->fdata = fftwf_alloc_real(coef->h * coef->w);
        if(!coef->fdata) { die("allocation error"); }
        unsigned blocks = (coef->h / 8) * (coef->w / 8);
        for(unsigned i = 0; i < blocks; i++) {
                for(unsigned j = 0; j < 64; j++) {
                        coef->fdata[i*64+j] = coef->data[i*64+j] * quant_table[j];
                }
                for(unsigned v = 0; v < 8; v++) {
                        for(unsigned u = 0; u < 8; u++) {
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
        for(unsigned i = 0; i < blocks * 64; i++) {
                coef->fdata[i] /= 16.;
        }
}
