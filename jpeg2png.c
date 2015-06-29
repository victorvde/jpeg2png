#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdnoreturn.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "gopt/gopt.h"

#include "jpeg2png.h"
#include "utils.h"
#include "jpeg.h"
#include "png.h"
#include "box.h"
#include "compute.h"
#include "logger.h"
#include "progressbar.h"
#include "fp_exceptions.h"

#define JPEG2PNG_VERSION "0.5"
static const float default_weight = 0.3;
static const float default_pweight = 0.001;
static const unsigned default_iterations = 50;

noreturn static void usage() {
        printf(
                "usage: jpeg2png picture.jpg ... [-o picture.png] ... [flags...]\n"
                "\n"
                "-o picture.png\n"
                "--output picture.png\n"
                "\tpicture.png is the file name of the output file\n"
                "\tthe output file will be overwritten if this flag is used\n"
                "\tmust be specified either zero times or once for every input file\n"
                "\tdefault value: original file name with the extension .png\n"
                "\n");
        printf(
                "-f\n"
                "--force\n"
                "\toverwrite output files even when not given explicit file names\n"
                "\n");
        printf(
                "-w weight[,weight_cb,weight_cr]\n"
                "--second-order-weight weight[,weight_cb,weight_cr]\n"
                "\tweight is a floating point number for TGV weight alpha_1\n"
                "\thigher values give smoother transitions with less staircasing\n"
                "\ta value of 1.0 means equivalent weight to the first order weight\n"
                "\ta value of 0.0 means plain Total Variation, and gives a speed boost\n"
                "\tweights for the chroma components always default to 0.\n"
                "\tdefault value: %g\n"
                "\n", default_weight);
        printf(
                "-p pweight[,pweight_cb,pweight_cr]\n"
                "--probability-weight pweight[,pweight_cb,pweight_cr]\n"
                "\tpweight is a floating point number for DCT coefficient distance weight\n"
                "\thigher values make the result more similar to the source JPEG\n"
                "\ta value of 1.0 means about equivalent weight to the first order weight\n"
                "\ta value of 0.0 means to ignore this and gives a speed boost\n"
                "\tweights for the chroma components default to the luma weight\n"
                "\tdefault value: %g\n"
                "\n", default_pweight);
        printf(
                "-i iterations[,iterations_cb,iterations_cr]\n"
                "--iterations iterations[,iterations_cb,iterations_cr]\n"
                "\titerations is an integer for the number of optimization steps\n"
                "\thigher values give better results but take more time\n"
                "\titerations for the chroma components default to the luma iterations\n"
                "\tdefault value: %d\n"
                "\n", default_iterations);
        printf(
                "-q\n"
                "--quiet\n"
                "\tdon't show the progress bar\n"
                "\n");
        printf(
                "-s\n"
                "--separate-components\n"
                "\tseparately optimize components\n"
                "\tthis is faster and makes multithreading more effective\n"
                "\thowever the edges of different components can be different\n"
                "\n");
        printf(
                "-t threads\n"
                "--threads threads\n"
#ifndef _OPENMP
                "\t*this version was compiled without support for threads*\n"
                "\n"
#endif
                "\tthreads is a positive integer for the maximum number of threads used\n"
                "\tequivalent to setting the environment variable OMP_NUM_THREADS\n"
                "\tdefault: number of CPUs\n"
                "\n");
        printf(
                "-1\n"
                "--16-bits-png\n"
                "\toutput PNG with 16 bits color depth instead of the usual 8 bits\n"
                "\tyou should use a high number of iterations when using this option\n"
                "\n");
        printf(
                "-c csv_log\n"
                "--csv_log csv_log\n"
                "\tcsv_log is a file name for the optimization log\n"
                "\tdefault: none\n"
                "\n");
        printf(
                "-h\n"
                "--help\n"
                "\tdisplay this help text and exit\n"
                "\n");
        printf(
                "-V\n"
                "--version\n"
                "\tdisplay version information and exit\n"
                );
        exit(EXIT_FAILURE);
}

void decode_file(const char* infile, const char *outfile, unsigned iterations[3], float weights[3], float pweights[3], unsigned png_bits, bool all_together, struct progressbar *pb, struct logger *plog) {
        FILE *in = fopen(infile, "rb");
        if(!in) { die_perror("could not open input file `%s`", infile); }
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

                unbox(coef->fdata, temp, coef->w, coef->h);

                free_real(coef->fdata);
                coef->fdata = temp;
        }

        if(all_together) {
                plog->channel = 3;
                compute(3, jpeg.coefs, plog, pb, weights[0], pweights, iterations[0]);
        } else {
                struct logger log = *plog;
                OPENMP(parallel for schedule(dynamic) firstprivate(log))
                for(unsigned i = 0; i < 3; i++) {
                        log.channel = i;
                        struct coef *coef = &jpeg.coefs[i];
                        compute(1, coef, &log, pb, weights[i], &pweights[i], iterations[i]);
                }
        }

        struct coef *coef = &jpeg.coefs[0];
        for(unsigned i = 0; i < coef->h * coef->w; i++) {
                coef->fdata[i] += 128.;
        }

        FILE *out = fopen(outfile, "wb");
        if(!out) { die_perror("could not open output file `%s`", outfile); }
        write_png(out, jpeg.w, jpeg.h, png_bits, &jpeg.coefs[0], &jpeg.coefs[1], &jpeg.coefs[2]);
        fclose(out);

        for(unsigned i = 0; i < 3; i++) {
                free_real(jpeg.coefs[i].fdata);
                free(jpeg.coefs[i].data);
        }
}


struct progressbar *main_progressbar;

int main(int argc, const char **argv) {
        enable_fp_exceptions();

        void *options = gopt_sort(&argc, argv, gopt_start(
                gopt_option('h', GOPT_NOARG, gopt_shorts('h','?'), gopt_longs("help")),
                gopt_option('V', GOPT_NOARG, gopt_shorts('V'), gopt_longs("version")),
                gopt_option('o', GOPT_ARG | GOPT_REPEAT, gopt_shorts('o'), gopt_longs("output")),
                gopt_option('f', GOPT_NOARG, gopt_shorts('f'), gopt_longs("force")),
                gopt_option('c', GOPT_ARG, gopt_shorts('c'), gopt_longs("csv-log")),
                gopt_option('t', GOPT_ARG, gopt_shorts('t'), gopt_longs("threads")),
                gopt_option('q', GOPT_NOARG, gopt_shorts('q'), gopt_longs("quiet")),
                gopt_option('s', GOPT_NOARG, gopt_shorts('s'), gopt_longs("separate-components")),
                gopt_option('1', GOPT_NOARG, gopt_shorts('1'), gopt_longs("16-bits-png")),
                gopt_option('i', GOPT_ARG, gopt_shorts('i'), gopt_longs("iterations")),
                gopt_option('p', GOPT_ARG, gopt_shorts('p'), gopt_longs("probability-weight")),
                gopt_option('w', GOPT_ARG, gopt_shorts('w'), gopt_longs("second-order-weight"))));
        if(gopt(options, 'V')) {
                printf("jpeg2png version "JPEG2PNG_VERSION" licensed GPLv3+\n");
                exit(EXIT_FAILURE);
        }
        if(argc < 2 || gopt(options, 'h')) {
                usage();
        }

        bool all_together = ! gopt(options, 's');

        const char *arg_string;
        float weights[3] = {default_weight, 0., 0.};
        if(gopt_arg(options, 'w', &arg_string)) {
                int n = sscanf(arg_string, "%f,%f,%f", &weights[0], &weights[1], &weights[2]);
                if(n == 3) {
                        if(all_together) {
                                die("different weights are only possible when using separated components");
                        }
                } else if(n ==1) {
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
                                die("different iteration counts are only possible when using separated components");
                        }
                } else if(n == 1) {
                        iterations[1] = iterations[0];
                        iterations[2] = iterations[0];
                } else {
                        die("invalid number of iterations");
                }
        }

        if(gopt_arg(options, 't', &arg_string)) {
#ifdef _OPENMP
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

        FILE *csv_log = NULL;
        if(gopt_arg(options, 'c', &arg_string)) {
                csv_log = fopen(arg_string, "wb");
                if(!csv_log) { die_perror("could not open csv log `%s`", csv_log); }
        }

        bool quiet = gopt(options, 'q');
        unsigned png_bits = gopt(options, '1') ? 16 : 8;
        bool force = gopt(options, 'f');

        struct logger log;
        logger_start(&log, csv_log);

        unsigned nin = argc - 1;
        unsigned nout = gopt(options, 'o');
        if(!(nout == 0 || nout == nin)) {
                die("must give output file names for all input files or none");
        }

        char **outfiles = malloc(sizeof(*outfiles) * nin);
        if(!outfiles) { die("could not allocate outfiles"); }
        if(nout) {
                gopt_args(options, 'o', (const char **)outfiles, nout);
        } else {
                for(unsigned i = 0; i < nin; i++) {
                        const char *infile = argv[1+i];
                        FILE *in = fopen(infile, "rb");
                        if(!in) { die("could not open input file `%s`", infile); }
                        fclose(in);

                        unsigned l = strlen(infile);
                        unsigned e = l;
                        if(l >= 5 && memcmp(".jpeg", &infile[l-5], 5) == 0) {
                                e = l-5;
                        } else if(l >= 4 && memcmp(".jpg", &infile[l-4], 4) == 0) {
                                e = l-4;
                        }
                        char *outfile = malloc(e + 4 + 1);
                        if(!outfile) { die("could not allocate outfile"); }
                        memcpy(outfile, infile, e);
                        memcpy(outfile+e, ".png", 5);

                        if(!nout && !force) {
                                // don't overwrite when not given -o or -f, racy
                                FILE *out = fopen(outfile, "rb");
                                if(out) { die("not overwriting output file `%s`", outfile); }
                        }

                        FILE *out = fopen(outfile, "wb");
                        if(!out) { die("could not open output file `%s`", outfile); }
                        fclose(out);
                        remove(outfile);

                        outfiles[i] = outfile;
               }
        }

        struct progressbar pb;
        if(!quiet) {
                if(all_together) {
                        progressbar_start(&pb, nin * iterations[0]);
                } else {
                        progressbar_start(&pb, nin * (iterations[0] + iterations[1] + iterations[2]));
                }
                main_progressbar = &pb;
        }

        OPENMP(parallel for schedule(dynamic) if(nin > 1) firstprivate(log))
        for(unsigned i = 0; i < nin; i++) {
                const char *infile = argv[1+i];
                const char *outfile = outfiles[i];
                log.filename = infile;

                decode_file(infile, outfile, iterations, weights, pweights, png_bits, all_together, quiet ? NULL : &pb, &log);
        }

        if(!nout) {
                for(unsigned i = 0; i < nin; i++) {
                        free(outfiles[i]);
                }
        }
        free(outfiles);

        if(!quiet) {
                progressbar_clear(&pb);
                main_progressbar = NULL;
        }

        gopt_free(options);

        if(csv_log) {
                fclose(csv_log);
        }
}
