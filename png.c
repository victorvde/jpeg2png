#include <stdio.h>
#include <stdlib.h>
#include <png.h>

#include "png.h"
#include "utils.h"

static float clamp(float x) {
        return CLAMP(x, 0., 255.);
}

void write_png(FILE *out, unsigned w, unsigned h, struct coef *y, struct coef *cb, struct coef *cr) {
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
        png_set_IHDR(png_ptr, info_ptr, w, h, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
        png_write_info(png_ptr, info_ptr);
        png_byte *image_data = calloc(sizeof(png_byte), h * w * 3);
        if(!image_data) { die("could not allocate image data");}

        for(unsigned i = 0; i < h; i++) {
                for(unsigned j = 0; j < w; j++) {
                        float yi = *p(y->fdata, i, j, y->w, y->h);
                        float cbi = *p(cb->fdata, i, j, cb->w, cb->h);
                        float cri = *p(cr->fdata, i, j, cr->w, cb->h);

                        image_data[(i*w+j)*3] = clamp(yi + 1.402 * cri);
                        image_data[(i*w+j)*3+1] = clamp(yi - 0.34414 * cbi - 0.71414 * cri);
                        image_data[(i*w+j)*3+2] = clamp(yi + 1.772 * cbi);
                }
        }

        png_byte **rows = malloc(sizeof(*rows) * h);
        if(!rows) { die("allocation failure"); }
        for(unsigned i = 0; i < h; i++) {
                rows[i] = &image_data[i * w * 3];
        }
        png_write_image(png_ptr, rows);
        free(rows);
        png_write_end(png_ptr, info_ptr);
        free(image_data);
        png_destroy_write_struct(&png_ptr, &info_ptr);
}
