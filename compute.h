#ifndef JPEG2PNG_COMPUTE_H
#define JPEG2PNG_COMPUTE_H

#include <stdint.h>
#include "logger.h"
#include "progressbar.h"

void compute(unsigned nchannel, struct coef coefs[nchannel], struct logger *log, struct progressbar *pb, float weight, float pweight[nchannel], unsigned iterations);

#endif
