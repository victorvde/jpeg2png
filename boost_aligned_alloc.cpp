#include <boost/align.hpp>
#include <boost/align/aligned_alloc.hpp>
#include "utils.h"

void *boost_aligned_alloc(size_t alignment, size_t size) {
	return boost::alignment::aligned_alloc(alignment, size);
}


void boost_aligned_free(void *p) {
	boost::alignment::aligned_free(p);
}

