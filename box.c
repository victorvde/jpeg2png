#include <assert.h>

#include "box.h"
#include "utils.h"

void unbox(float *restrict in, float *restrict out, unsigned w, unsigned h) {
        assert((w & 7) == 0);
        assert((h & 7) == 0);
        for(unsigned block_y = 0; block_y < h / 8; block_y++) {
                for(unsigned block_x = 0; block_x < w / 8; block_x++) {
                        for(unsigned in_y = 0; in_y < 8; in_y++) {
                                for(unsigned in_x = 0; in_x < 8; in_x++) {
                                        *p(out, block_x * 8 + in_x, block_y * 8 + in_y, w, h) = *in++;
                                }
                        }
                }
        }
}

void box(float *restrict in, float *restrict out, unsigned w, unsigned h) {
        assert((w & 7) == 0);
        assert((h & 7) == 0);
        for(unsigned block_y = 0; block_y < h / 8; block_y++) {
                for(unsigned block_x = 0; block_x < w / 8; block_x++) {
                        for(unsigned in_y = 0; in_y < 8; in_y++) {
                                for(unsigned in_x = 0; in_x < 8; in_x++) {
                                        *out++ = *p(in, block_x * 8 + in_x, block_y * 8 + in_y, w, h);
                                }
                        }
                }
        }
}
