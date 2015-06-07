#include <stdio.h>
#include "progressbar.h"

static const unsigned progressbar_width = 70;

void progressbar_start(struct progressbar *pb, unsigned max) {
        pb->max = max;
        progressbar_set(pb, 0);
}

void progressbar_set(struct progressbar *pb, unsigned current) {
        unsigned old_to_print = progressbar_width * pb->current / pb->max;
        unsigned old_percentage = 100 * pb->current / pb->max;
        pb->current = current;
        unsigned to_print = progressbar_width * pb->current / pb->max;
        unsigned percentage = 100 * pb->current / pb->max;
        if(old_to_print == to_print && old_percentage == percentage) {
                return;
        }
        printf("\r[");
        for(unsigned i = 0; i < to_print; i++) {
                printf("#");
        }
        for(unsigned i = 0; i < progressbar_width - to_print; i++) {
                printf(" ");
        }
        printf("] %3d%%", percentage);
        fflush(stdout);
}

void progressbar_add(struct progressbar *pb, unsigned n) {
        progressbar_set(pb, pb->current + n);
}

void progressbar_inc(struct progressbar *pb) {
        progressbar_add(pb, 1);
}

void progressbar_done(struct progressbar *pb) {
        (void)pb;
        printf("\r");
        for(unsigned i = 0; i < progressbar_width + 7; i++) {
                printf(" ");
        }
        fflush(stdout);
        printf("\r");
        fflush(stdout);
}
