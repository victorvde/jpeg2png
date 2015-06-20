#include <stdio.h>

#include "logger.h"
#include "utils.h"

void logger_start(struct logger *log, FILE *csv_log) {
        log->f = csv_log;
        log->filename = "";
        log->channel = 0;
        log->iteration = 0;
        if(log->f) {
                if(fprintf(log->f, "filename,channel,iteration,objective,prob_dist,tv,tv2\n") < 0) {
                        die_perror("could not write to csv log");
                }
        }
}

void logger_log(struct logger *log, double objective, double prob_dist, double tv, double tv2) {
        if(log->f) {
#ifdef USE_OPENMP
        #pragma omp critical(write_log)
#endif
                if(fprintf(log->f, "%s,%u,%u,%f,%f,%f,%f\n", log->filename, log->channel, log->iteration, objective, prob_dist, tv, tv2) < 0)
                {
                        die_perror("could not write to csv log");
                }
        }
}
