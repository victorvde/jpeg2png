#include <stdio.h>

#include "logger.h"
#include "utils.h"

void logger_start(struct logger *log, FILE *csv_log) {
        log->f = csv_log;
        log->channel = -1;
        if(log->f) {
                if(fprintf(log->f, "channel,iteration,objective,tv,tv2\n") < 0) {
                        die_perror("could not write to csv log");
                }
        }
}

void logger_log(struct logger *log, double objective, double tv, double tv2) {
        if(log->f) {
#ifdef USE_OPENMP
        #pragma omp critical(write_log)
#endif
                if(fprintf(log->f, "%d,%d,%f,%f,%f\n", log->channel, log->iteration, objective, tv, tv2) < 0)
                {
                        die_perror("could not write to csv log");
                }
        }
}
