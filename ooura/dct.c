/*
Source: http://www.kurims.kyoto-u.ac.jp/~ooura/fft.html
fft2d.zip (2006/12/28) file shrtdct.c

Copyright Takuya OOURA, 1996-2001

You may use, copy, modify and distribute this code for any purpose
(include commercial use) and without fee. Please refer to this package
when you modify this code.

Modifications:
double -> float
double indirection removed
added named functions instead of sign
*/

// Cn_kR = sqrt(2.0/n) * cos(pi/2*k/n)
// Cn_kI = sqrt(2.0/n) * sin(pi/2*k/n)
// Wn_kR = cos(pi/2*k/n)
// Wn_kI = sin(pi/2*k/n)
#define C8_1R   0.49039264020161522456
#define C8_1I   0.09754516100806413392
#define C8_2R   0.46193976625564337806
#define C8_2I   0.19134171618254488586
#define C8_3R   0.41573480615127261854
#define C8_3I   0.27778511650980111237
#define C8_4R   0.35355339059327376220
#define W8_4R   0.70710678118654752440

// Normalized 8x8 IDCT
void idct8x8s(float a[64]) {
        float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
        float xr, xi;

        for (unsigned j = 0; j <= 7; j++) {
                x1r = C8_1R * a[1*8+j] + C8_1I * a[7*8+j];
                x1i = C8_1R * a[7*8+j] - C8_1I * a[1*8+j];
                x3r = C8_3R * a[3*8+j] + C8_3I * a[5*8+j];
                x3i = C8_3R * a[5*8+j] - C8_3I * a[3*8+j];
                xr = x1r - x3r;
                xi = x1i + x3i;
                x1r += x3r;
                x3i -= x1i;
                x1i = W8_4R * (xr + xi);
                x3r = W8_4R * (xr - xi);
                xr = C8_2R * a[2*8+j] + C8_2I * a[6*8+j];
                xi = C8_2R * a[6*8+j] - C8_2I * a[2*8+j];
                x0r = C8_4R * (a[0*8+j] + a[4*8+j]);
                x0i = C8_4R * (a[0*8+j] - a[4*8+j]);
                x2r = x0r - xr;
                x2i = x0i - xi;
                x0r += xr;
                x0i += xi;
                a[0*8+j] = x0r + x1r;
                a[7*8+j] = x0r - x1r;
                a[2*8+j] = x0i + x1i;
                a[5*8+j] = x0i - x1i;
                a[4*8+j] = x2r - x3i;
                a[3*8+j] = x2r + x3i;
                a[6*8+j] = x2i - x3r;
                a[1*8+j] = x2i + x3r;        
        }
        for (unsigned j = 0; j <= 7; j++) {
                x1r = C8_1R * a[j*8+1] + C8_1I * a[j*8+7];
                x1i = C8_1R * a[j*8+7] - C8_1I * a[j*8+1];
                x3r = C8_3R * a[j*8+3] + C8_3I * a[j*8+5];
                x3i = C8_3R * a[j*8+5] - C8_3I * a[j*8+3];
                xr = x1r - x3r;
                xi = x1i + x3i;
                x1r += x3r;
                x3i -= x1i;
                x1i = W8_4R * (xr + xi);
                x3r = W8_4R * (xr - xi);
                xr = C8_2R * a[j*8+2] + C8_2I * a[j*8+6];
                xi = C8_2R * a[j*8+6] - C8_2I * a[j*8+2];
                x0r = C8_4R * (a[j*8+0] + a[j*8+4]);
                x0i = C8_4R * (a[j*8+0] - a[j*8+4]);
                x2r = x0r - xr;
                x2i = x0i - xi;
                x0r += xr;
                x0i += xi;
                a[j*8+0] = x0r + x1r;
                a[j*8+7] = x0r - x1r;
                a[j*8+2] = x0i + x1i;
                a[j*8+5] = x0i - x1i;
                a[j*8+4] = x2r - x3i;
                a[j*8+3] = x2r + x3i;
                a[j*8+6] = x2i - x3r;
                a[j*8+1] = x2i + x3r;
        }
}

// Normalized 8x8 DCT
void dct8x8s(float a[64]) {
        float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;
        float xr, xi;

        for (unsigned j = 0; j <= 7; j++) {
                x0r = a[0*8+j] + a[7*8+j];
                x1r = a[0*8+j] - a[7*8+j];
                x0i = a[2*8+j] + a[5*8+j];
                x1i = a[2*8+j] - a[5*8+j];
                x2r = a[4*8+j] + a[3*8+j];
                x3r = a[4*8+j] - a[3*8+j];
                x2i = a[6*8+j] + a[1*8+j];
                x3i = a[6*8+j] - a[1*8+j];
                xr = x0r + x2r;
                xi = x0i + x2i;
                a[0*8+j] = C8_4R * (xr + xi);
                a[4*8+j] = C8_4R * (xr - xi);
                xr = x0r - x2r;
                xi = x0i - x2i;
                a[2*8+j] = C8_2R * xr - C8_2I * xi;
                a[6*8+j] = C8_2R * xi + C8_2I * xr;
                xr = W8_4R * (x1i - x3i);
                x1i = W8_4R * (x1i + x3i);
                x3i = x1i - x3r;
                x1i += x3r;
                x3r = x1r - xr;
                x1r += xr;
                a[1*8+j] = C8_1R * x1r - C8_1I * x1i;
                a[7*8+j] = C8_1R * x1i + C8_1I * x1r;
                a[3*8+j] = C8_3R * x3r - C8_3I * x3i;
                a[5*8+j] = C8_3R * x3i + C8_3I * x3r;
        }
        for (unsigned j = 0; j <= 7; j++) {
                x0r = a[j*8+0] + a[j*8+7];
                x1r = a[j*8+0] - a[j*8+7];
                x0i = a[j*8+2] + a[j*8+5];
                x1i = a[j*8+2] - a[j*8+5];
                x2r = a[j*8+4] + a[j*8+3];
                x3r = a[j*8+4] - a[j*8+3];
                x2i = a[j*8+6] + a[j*8+1];
                x3i = a[j*8+6] - a[j*8+1];
                xr = x0r + x2r;
                xi = x0i + x2i;
                a[j*8+0] = C8_4R * (xr + xi);
                a[j*8+4] = C8_4R * (xr - xi);
                xr = x0r - x2r;
                xi = x0i - x2i;
                a[j*8+2] = C8_2R * xr - C8_2I * xi;
                a[j*8+6] = C8_2R * xi + C8_2I * xr;
                xr = W8_4R * (x1i - x3i);
                x1i = W8_4R * (x1i + x3i);
                x3i = x1i - x3r;
                x1i += x3r;
                x3r = x1r - xr;
                x1r += xr;
                a[j*8+1] = C8_1R * x1r - C8_1I * x1i;
                a[j*8+7] = C8_1R * x1i + C8_1I * x1r;
                a[j*8+3] = C8_3R * x3r - C8_3I * x3i;
                a[j*8+5] = C8_3R * x3i + C8_3I * x3r;
        }
}
