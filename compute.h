#ifndef JPEG2PNG_COMPUTE_H
#define JPEG2PNG_COMPUTE_H

#include <stdint.h>
#include "logger.h"

void compute(struct coef *coef, struct logger *log, uint16_t quant_table[64], float weight, unsigned iterations);

#endif
