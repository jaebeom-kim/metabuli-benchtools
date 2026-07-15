// Replacement for mmseqs' FastSort.h. Upstream optionally uses ips4o; here we
// always fall back to std::sort, which is all StringBlock::compact() needs.
#ifndef COMPAT_FASTSORT_H
#define COMPAT_FASTSORT_H

#include <algorithm>

#define SORT_SERIAL std::sort
#define SORT_PARALLEL std::sort

#endif // COMPAT_FASTSORT_H
