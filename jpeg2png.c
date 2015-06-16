#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdnoreturn.h>
#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "gopt/gopt.h"

#include "jpeg2png.h"
#include "utils.h"
#include "jpeg.h"
#include "png.h"
#include "box.h"
#include "upsample.h"
#include "compute.h"
#include "logger.h"
#include "progressbar.h"

#define JPEG2PNG_VERSION "0.2"
static const float default_weight = 0.3;
static const float default_pweight = 0.001;
static const unsigned default_iterations = 50;

noreturn static void usage() {
        printf(
                "usage: jpeg2png in.jpg out.png [flags...]\n"
                "\n"
                "-w weight[,weight_cb,weight_cr]\n"
                "--second-order-weight weight[,weight_cb,weight_cr]\n"
                "\tweight is a floating point number for TVG weight alpha_1\n"
                "\thigher values give smoother transitions with less staircasing\n"
                "\ta value of 1.0 means equivalent weight to the first order weight\n"
                "\ta value of 0.0 means plain Total Variation, and gives a speed boost\n"
                "\tweights for the chroma components always default to 0.\n"
                "\tdefault value: %g\n"
                "\n"
                "-p pweight[,pweight_cb,pweight_cr]\n"
                "--probability-weight pweight[,pweight_cb,pweight_cr]\n"
                "\tpweight is a floating point number for DCT coefficient distance weight\n"
                "\thigher values make the result more similar to the source JPEG\n"
                "\ta value of 1.0 means about equivalent weight to the first order weight\n"
                "\ta value of 0.0 means to ignore this and gives a speed boost\n"
                "\tweights for the chroma components default to the luma weight\n"
                "\tdefault value: %g\n"
                "\n"
                "-i iterations[,iterations_cb,iterations_cr]\n"
                "--iterations iterations[,iterations_cb,iterations_cr]\n"
                "\titerations is an integer for the number of optimization steps\n"
                "\thigher values give better results but take more time\n"
                "\titerations for the chroma components default to the luma iterations\n"
                "\tdefault value: %d\n"
                "\n"
                "-q\n"
                "--quiet\n"
                "\tdon't show the progress bar\n"
                "\n"
                "-s\n"
                "--seperate-components\n"
                "\tseperately optimize components\n"
                "\tthis is faster and makes multithreading more effective\n"
                "\thowever the edges of different components can be different\n"
                "\n"
                "-t threads\n"
                "--threads threads\n"
#ifndef USE_OPENMP
                "\t*this version was compiled without support for threads*\n"
                "\n"
#endif
                "\tthreads is a positive integer for the maximum number of threads used\n"
                "\tequivalent to setting the environment variable OMP_NUM_THREADS\n"
                "\tdefault: number of CPUs\n"
                "\n"
                "-1\n"
                "--16-bits-png\n"
                "\toutput PNG with 16 bits color depth instead of the usual 8 bits\n"
                "\tyou should use a high number of iterations when using this option\n"
                "\n"
                "-c csv_log\n"
                "--csv_log csv_log\n"
                "\tcsv_log is a file name for the optimization log\n"
                "\tdefault: none\n"
                "\n"
                "-h\n"
                "--help\n"
                "\tdisplay this help text and exit\n"
                "\n"
                "-V\n"
                "--version\n"
                "\tdisplay version information and exit\n"
                , default_weight, default_pweight, default_iterations);
        exit(EXIT_FAILURE);
}

int main(int argc, const char **argv) {
        void *options = gopt_sort(&argc, argv, gopt_start(
                gopt_option('h', GOPT_NOARG, gopt_shorts('h','?'), gopt_longs("help")),
                gopt_option('V', GOPT_NOARG, gopt_shorts('V'), gopt_longs("version")),
                gopt_option('c', GOPT_ARG, gopt_shorts('c'), gopt_longs("csv-log")),
                gopt_option('t', GOPT_ARG, gopt_shorts('t'), gopt_longs("threads")),
                gopt_option('q', GOPT_NOARG, gopt_shorts('q'), gopt_longs("quiet")),
                gopt_option('s', GOPT_NOARG, gopt_shorts('s'), gopt_longs("seperate-components")),
                gopt_option('1', GOPT_NOARG, gopt_shorts('1'), gopt_longs("16-bits-png")),
                gopt_option('i', GOPT_ARG, gopt_shorts('i'), gopt_longs("iterations")),
                gopt_option('p', GOPT_ARG, gopt_shorts('p'), gopt_longs("probability-weight")),
                gopt_option('w', GOPT_ARG, gopt_shorts('w'), gopt_longs("second-order-weight"))));
        if(gopt(options, 'V')) {
                printf("jpeg2png version "JPEG2PNG_VERSION" licensed GPLv3+\n");
                exit(EXIT_FAILURE);
        }
        if(argc != 3 || gopt(options, 'h')) {
                usage();
        }

        bool all_together = ! gopt(options, 's');

        const char *arg_string;
        float weights[3] = {default_weight, 0., 0.};
        if(gopt_arg(options, 'w', &arg_string)) {
                int n = sscanf(arg_string, "%f,%f,%f", &weights[0], &weights[1], &weights[2]);
                if(n == 3 || n == 1) {
                        // ok
                } else {
                        die("invalid weight");
                }
        }
        float pweights[3] = {default_pweight, default_pweight, default_pweight};
        if(gopt_arg(options, 'p', &arg_string)) {
                int n = sscanf(arg_string, "%f,%f,%f", &pweights[0], &pweights[1], &pweights[2]);
                if(n == 3) {
                        // ok
                } else if(n == 1) {
                        pweights[1] = pweights[0];
                        pweights[2] = pweights[0];
                } else {
                        die("invalid probability weight");
                }
        }
        unsigned iterations[3] = {default_iterations, default_iterations, default_iterations};
        if(gopt_arg(options, 'i', &arg_string)) {
                int n = sscanf(arg_string, "%u,%u,%u", &iterations[0], &iterations[1], &iterations[2]);
                if(n == 3) {
                        if(all_together) {
                                die("different iteration counts are only possible when using seperated components");
                        }
                } else if(n == 1) {
                        iterations[1] = iterations[0];
                        iterations[2] = iterations[0];
                } else {
                        die("invalid number of iterations");
                }
        }

        if(gopt_arg(options, 't', &arg_string)) {
#ifdef USE_OPENMP
                unsigned threads;
                int n = sscanf(arg_string, "%u", &threads);
                if(n != 1 || threads == 0) {
                        die("invalid number of threads");
                }
                omp_set_num_threads(threads);
#else
                die("this version is compiled without support for threads");
#endif
        }

        FILE *in = fopen(argv[1], "rb");
        if(!in) { die_perror("could not open input file `%s`", argv[1]); }
        FILE *out = fopen(argv[2], "wb");
        if(!out) { die_perror("could not open output file `%s`", argv[2]); }

        FILE *csv_log = NULL;
        if(gopt_arg(options, 'c', &arg_string)) {
                csv_log = fopen(arg_string, "wb");
                if(!csv_log) { die_perror("could not open csv log `%s`", csv_log); }
        }

        bool quiet = gopt(options, 'q');
        unsigned png_bits = gopt(options, '1') ? 16 : 8;

        gopt_free(options);

        struct jpeg jpeg;
        read_jpeg(in, &jpeg);
        fclose(in);

        for(unsigned c = 0; c < 3; c++) {
                struct coef *coef = &jpeg.coefs[c];
                decode_coefficients(coef);
        }

        for(unsigned i = 0; i < 3; i++) {
                struct coef *coef = &jpeg.coefs[i];
                float *temp = alloc_real(coef->h * coef->w);
                if(!temp) { die("allocation error"); }

                unbox(coef->fdata, temp, coef->w, coef->h);

                free_real(coef->fdata);
                coef->fdata = temp;
        }

        struct logger log;
        logger_start(&log, csv_log);
        struct progressbar pb;
        if(all_together) {
                if(!quiet) {
                        progressbar_start(&pb, iterations[0]);
                }
                log.channel = 3;
                compute(3, jpeg.coefs, &log, quiet ? NULL : &pb, weights, pweights, iterations[0]);
        } else {
                if(!quiet) {
                        progressbar_start(&pb, iterations[0] + iterations[1] + iterations[2]);
                }
                #ifdef USE_OPENMP
                #pragma omp parallel for schedule(dynamic) firstprivate(log)
                #endif
                for(unsigned i = 0; i < 3; i++) {
                        log.channel = i;
                        struct coef *coef = &jpeg.coefs[i];
                        compute(1, coef, &log, quiet ? NULL : &pb, &weights[i], &pweights[i], iterations[i]);
                }
        }
        if(!quiet) {
                progressbar_done(&pb);
        }
        if(csv_log) {
                fclose(csv_log);
        }

        struct coef *coef = &jpeg.coefs[0];
        for(unsigned i = 0; i < coef->h * coef->w; i++) {
                coef->fdata[i] += 128.;
        }

        for(unsigned i = 0; i < 3; i++) {
                upsample(&jpeg.coefs[i], jpeg.w, jpeg.h);
        }

        write_png(out, jpeg.w, jpeg.h, png_bits, &jpeg.coefs[0], &jpeg.coefs[1], &jpeg.coefs[2]);
        fclose(out);

        for(unsigned i = 0; i < 3; i++) {
                free_real(jpeg.coefs[i].fdata);
                free(jpeg.coefs[i].data);
        }
}
