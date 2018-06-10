#ifndef JPEG2PNG_PROGRESSBAR_H
#define JPEG2PNG_PROGRESSBAR_H

struct progressbar {
        unsigned current;
        unsigned max;
};

// initialize progressbar with maximum value
void progressbar_start(struct progressbar *pb, unsigned max);
// set current value
void progressbar_set(struct progressbar *pb, unsigned current);
// add to current value
void progressbar_add(struct progressbar *pb, unsigned n);
// add one to current value
void progressbar_inc(struct progressbar *pb);
// clear line of progress bar
void progressbar_clear(struct progressbar *pb);

#endif
