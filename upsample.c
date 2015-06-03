#include <assert.h>

#include <fftw3.h>

#include "jpeg2png.h"

void upsample(struct coef *coef, int w, int h) {
        if(coef->h < h) {
                assert(coef->h / 8 == ((h + 7) / 8 + 1) / 2);
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
        if(coef->w < w) {
                assert(coef->w / 8 == ((w + 7) / 8 + 1) / 2);
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
