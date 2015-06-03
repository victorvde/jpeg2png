#ifndef JPEG2PNG_PNG_H
#define JPEG2PNG_PNG_H

#include <stdio.h>
#include "jpeg2png.h"

void write_png(FILE *out, int w, int h, struct coef *y, struct coef *cb, struct coef *cr);

#endif
