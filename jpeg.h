#ifndef JPEG2PNG_JPEG_H
#define JPEG2PNG_JPEG_H

#include <stdio.h>

#include "jpeg2png.h"

struct jpeg {
        int h;
        int w;
        uint16_t quant_table[3][64];
        struct coef coefs[3];
};

void read_jpeg(FILE *in, struct jpeg *jpeg);
#endif
