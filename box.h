#ifndef JPEG2PNG_BOX_H
#define JPEG2PNG_BOX_H

void unbox(float *restrict in, float *restrict out, unsigned w, unsigned h);
void box(float *restrict in, float *restrict out, unsigned w, unsigned h);

#endif
