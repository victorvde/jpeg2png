#ifndef JPEG2PNG_COMPUTE_H
#define JPEG2PNG_COMPUTE_H

#include <stdint.h>
#include "logger.h"
#include "progressbar.h"

void compute(struct coef *coef, struct logger *log, struct progressbar *pb, uint16_t quant_table[64], float weight, float pweight, unsigned iterations);

#endif
