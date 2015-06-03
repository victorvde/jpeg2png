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
#include "box.h"
#include "upsample.h"
#include "compute.h"

int main(int argc, char *argv[]) {
        if(argc != 3) {
                die("usage: jpeg2png in.jpg out.png");
        }

        FILE *in = fopen(argv[1], "rb");
        if(!in) { die(NULL); }
        FILE *out = fopen(argv[2], "wb");
        if(!out) { die(NULL); }

        struct jpeg jpeg;
        read_jpeg(in, &jpeg);
        fclose(in);

        for(int c = 0; c < 3; c++) {
                struct coef *coef = &jpeg.coefs[c];
                decode_coefficients(coef, jpeg.quant_table[c]);
        }

        for(int i = 0; i < 3; i++) {
                struct coef *coef = &jpeg.coefs[i];
                float *temp = fftwf_alloc_real(coef->h * coef->w);
                if(!temp) { die("allocation error"); }

                unbox(coef->fdata, temp, coef->w, coef->h);

                fftwf_free(coef->fdata);
                coef->fdata = temp;
        }

        START_TIMER(computing);
        for(int i = 0; i < 3; i++) {
                START_TIMER(compute_1);
                struct coef *coef = &jpeg.coefs[i];
                uint16_t *quant_table = jpeg.quant_table[i];
                compute(coef, quant_table);
                STOP_TIMER(compute_1);
        }
        STOP_TIMER(computing);

        struct coef *coef = &jpeg.coefs[0];
        for(int i = 0; i < coef->h * coef->w; i++) {
                coef->fdata[i] += 128.;
        }

        START_TIMER(upsampling);
        for(int i = 0; i < 3; i++) {
                upsample(&jpeg.coefs[i], jpeg.w, jpeg.h);
        }
        STOP_TIMER(upsampling);

        write_png(out, jpeg.w, jpeg.h, &jpeg.coefs[0], &jpeg.coefs[1], &jpeg.coefs[2]);
        fclose(out);

        for(int i = 0; i < 3; i++) {
                fftwf_free(jpeg.coefs[i].fdata);
                free(jpeg.coefs[i].data);
        }
}
