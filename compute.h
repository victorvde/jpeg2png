#ifndef JPEG2PNG_COMPUTE_H
#define JPEG2PNG_COMPUTE_H

#include <stdint.h>
#include "logger.h"
#include "progressbar.h"

void compute(unsigned ncoef, struct coef coefs[ncoef], struct logger *log, struct progressbar *pb, float weight[ncoef], float pweight[ncoef], unsigned iterations);

#endif
