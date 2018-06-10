#ifndef JPEG2PNG_JPEG2PNG_H
#define JPEG2PNG_JPEG2PNG_H

#include <stdint.h>
#include "progressbar.h"

struct coef {
        unsigned h;
        unsigned w;
        // vertical subsampling factor
        unsigned h_samp;
        // horizontal subsampling factor
        unsigned w_samp;
        // DCT coefficients
        int16_t *data;
        // image data
        float *fdata;
        // quantization table
        uint16_t quant_table[64];
};

extern struct progressbar *main_progressbar;

#endif
