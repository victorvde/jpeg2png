#ifndef JPEG2PNG_LOGGER_H
#define JPEG2PNG_LOGGER_H

#include <stdio.h>

struct logger {
        FILE *f;
        unsigned channel;
        unsigned iteration;
};

void logger_start(struct logger *log, FILE *csv_log);
void logger_log(struct logger *log, float objective, float tv, float tv2);

#endif
