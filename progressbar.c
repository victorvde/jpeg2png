#include <stdio.h>
#include "progressbar.h"

// see progressbar.h

static const unsigned progressbar_width = 70;

static unsigned get_to_print(unsigned current, unsigned max) {
        return progressbar_width * current / max;
}

static unsigned get_percentage(unsigned current, unsigned max) {
        return 100 * current / max;
}

static void progressbar_show(struct progressbar *pb) {
        unsigned to_print = get_to_print(pb->current, pb->max);
        unsigned percentage = get_percentage(pb->current, pb->max);

        printf("\r[");
        for(unsigned i = 0; i < to_print; i++) {
                printf("#");
        }
        for(unsigned i = 0; i < progressbar_width - to_print; i++) {
                printf(" ");
        }
        printf("] %3u%%", percentage);
        fflush(stdout);
}

void progressbar_start(struct progressbar *pb, unsigned max) {
        pb->max = max;
        pb->current = 0;
        progressbar_show(pb);
}

void progressbar_set(struct progressbar *pb, unsigned current) {
        unsigned old_to_print = get_to_print(pb->current, pb->max);
        unsigned old_percentage = get_percentage(pb->current, pb->max);
        unsigned to_print = get_to_print(current, pb->max);
        unsigned percentage = get_percentage(current, pb->max);
        pb->current = current;
        if(old_to_print == to_print && old_percentage == percentage) {
                return;
        }
        progressbar_show(pb);
}

void progressbar_add(struct progressbar *pb, unsigned n) {
        progressbar_set(pb, pb->current + n);
}

void progressbar_inc(struct progressbar *pb) {
        progressbar_add(pb, 1);
}

void progressbar_clear(struct progressbar *pb) {
        (void)pb;
        printf("\r");
        for(unsigned i = 0; i < progressbar_width + 7; i++) {
                printf(" ");
        }
        fflush(stdout);
        printf("\r");
        fflush(stdout);
}
