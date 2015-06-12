#include <stdio.h>
#include <stdlib.h>
#include <png.h>

#include "png.h"
#include "utils.h"

static float clamp(float x) {
        return CLAMP(x, 0., 255.);
}

void write_png(FILE *out, unsigned w, unsigned h, unsigned bits, struct coef *y, struct coef *cb, struct coef *cr) {
        ASSUME(bits == 8 || bits == 16);
        png_struct *png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
        if(!png_ptr) { die("could not initialize PNG write struct"); }
        png_info *info_ptr = png_create_info_struct(png_ptr);
        if (!info_ptr) { die("could not initialize PNG info struct"); }
        if (setjmp(png_jmpbuf(png_ptr)))
        {
                /* If we get here, we had a problem writing the file */
                die("problem while writing PNG file");
        }
        png_init_io(png_ptr, out);
        png_set_IHDR(png_ptr, info_ptr, w, h, bits, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_write_info(png_ptr, info_ptr);
        unsigned depth = bits / 8;
        png_byte *image_data = calloc(sizeof(png_byte), h * w * 3 * depth);
        if(!image_data) { die("could not allocate image data");}

        for(unsigned i = 0; i < h; i++) {
                for(unsigned j = 0; j < w; j++) {
                        float yi = *p(y->fdata, j, i, y->w, y->h);
                        float cbi = *p(cb->fdata, j, i, cb->w, cb->h);
                        float cri = *p(cr->fdata, j, i, cr->w, cr->h);

                        float bitfactor = (1 << bits) / 256.;
                        unsigned r = clamp(yi + 1.402 * cri) * bitfactor;
                        unsigned g = clamp(yi - 0.34414 * cbi - 0.71414 * cri) * bitfactor;
                        unsigned b = clamp(yi + 1.772 * cbi) * bitfactor;

                        png_byte *here = &image_data[(i*w+j)*3*depth];
                        if(bits == 8) {
                                here[0] = r & 0xFF;
                                here[1] = g & 0xFF;
                                here[2] = b & 0xFF;
                        } else {
                                here[0] = (r >> 8) & 0xFF;
                                here[1] = r & 0xFF;
                                here[2] = (g >> 8) & 0xFF;
                                here[3] = g & 0xFF;
                                here[4] = (b >> 8) & 0xFF;
                                here[5] = b & 0xFF;
                        }
                }
        }

        png_byte **rows = malloc(sizeof(*rows) * h);
        if(!rows) { die("allocation failure"); }
        for(unsigned i = 0; i < h; i++) {
                rows[i] = &image_data[i * w * 3 * depth];
        }
        png_write_image(png_ptr, rows);
        free(rows);
        png_write_end(png_ptr, info_ptr);
        free(image_data);
        png_destroy_write_struct(&png_ptr, &info_ptr);
}
