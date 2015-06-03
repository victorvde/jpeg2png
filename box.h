#ifndef JPEG2PNG_BOX_H
#define JPEG2PNG_BOX_H

void unbox(float *restrict in, float *restrict out, int w, int h);
void box(float *restrict in, float *restrict out, int w, int h);

#endif
