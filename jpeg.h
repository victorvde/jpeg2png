#ifndef JPEG2PNG_JPEG_H
#define JPEG2PNG_JPEG_H

#include <stdio.h>
#include <stdint.h>

#include "jpeg2png.h"

struct jpeg {
        unsigned c;
        unsigned h;
        unsigned w;
        struct coef coefs[4];
};

void read_jpeg(FILE *in, struct jpeg *jpeg);
void decode_coefficients(struct coef *coef);
#endif
