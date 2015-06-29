#ifndef JPEG2PNG_PROGRESSBAR_H
#define JPEG2PNG_PROGRESSBAR_H

struct progressbar {
        unsigned current;
        unsigned max;
};

void progressbar_start(struct progressbar *pb, unsigned max);
void progressbar_set(struct progressbar *pb, unsigned current);
void progressbar_add(struct progressbar *pb, unsigned n);
void progressbar_inc(struct progressbar *pb);
void progressbar_clear(struct progressbar *pb);

#endif
