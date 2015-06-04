#include <assert.h>

#include <fftw3.h>

#include "jpeg2png.h"
#include "utils.h"

static void upsample_y(struct coef *coef, int jpeg_w, int jpeg_h) {
        (void) jpeg_w;

        int w = coef->w;
        int h = coef->h;

        int jpeg_blocks_h = (jpeg_h + 7) / 8;
        assert(h / 8 == (jpeg_blocks_h + 1) / 2);
        (void) jpeg_blocks_h;

        int new_h = h * 2;
        float *new = fftwf_alloc_real(w * new_h);
        for(int y = 0; y < new_h; y++) {
                for(int x = 0; x < w; x++) {
                        if(y == 0) {
                                *p(new, x, y, w, new_h) = *p(coef->fdata, x, y / 2, w, h);
                        } else if (y == new_h-1) {
                                *p(new, x, y, w, new_h) = *p(coef->fdata, x, y / 2 - 1, w, h);
                        } else {
                                int other_pixel = y % 2 == 0 ? -1 : 1;
                                *p(new, x, y, w, new_h) =
                                        0.75 * *p(coef->fdata, x, y / 2, w, h) +
                                        0.25 * *p(coef->fdata, x, y / 2 + other_pixel, w, h);
                        }
                }
        }
        fftwf_free(coef->fdata);
        coef->fdata = new;
        coef->h = new_h;
}

static void upsample_x(struct coef *coef, int jpeg_w, int jpeg_h) {
        (void) jpeg_h;

        int w = coef->w;
        int h = coef->h;

        int jpeg_blocks_w = (jpeg_w + 7) / 8;
        assert(w / 8 == (jpeg_blocks_w + 1) / 2);
        (void) jpeg_blocks_w;

        int new_w = w * 2;
        float *new = fftwf_alloc_real(new_w * h);
        for(int y = 0; y < h; y++) {
                for(int x = 0; x < new_w; x++) {
                        if(x == 0) {
                                *p(new, x, y, new_w, h) = *p(coef->fdata, x / 2, y, w, h);
                        } else if (x == new_w-1) {
                                *p(new, x, y, new_w, h) = *p(coef->fdata, x / 2 - 1, y, w, h);
                        } else {
                                int other_pixel = x % 2 == 0 ? -1 : 1;
                                *p(new, x, y, new_w, h) =
                                        0.75 * *p(coef->fdata, x / 2, y, w, h) +
                                        0.25 * *p(coef->fdata, x / 2 + other_pixel, y, w, h);
                        }
                }
        }
        fftwf_free(coef->fdata);
        coef->fdata = new;
        coef->w = new_w;
}

void upsample(struct coef *coef, int jpeg_w, int jpeg_h) {
        if(coef->h < jpeg_h) {
                upsample_y(coef, jpeg_w, jpeg_h);
        }
        if(coef->w < jpeg_w) {
                upsample_x(coef, jpeg_w, jpeg_h);
        }
}
