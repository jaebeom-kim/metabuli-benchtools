// Minimal subset of mmseqs' MathUtil.h. Only flog2 and ipow are needed by the
// ported taxonomy code (sparse-table sizing for the LCA range-minimum query).
// The exact fast-log2 approximation from mmseqs is preserved so the sparse
// table is sized and indexed identically to upstream Metabuli.
#ifndef COMPAT_MATHUTIL_H
#define COMPAT_MATHUTIL_H

#if defined(__has_attribute)
#  if __has_attribute(__may_alias__)
#    define MAY_ALIAS(x) x __attribute__((__may_alias__))
#  else
#    define MAY_ALIAS(x) x
#  endif
#else
#  define MAY_ALIAS(x) x
#endif

class MathUtil {
public:
    template<typename T>
    static inline T ipow(int base, int exponent) {
        T res = 1;
        for (int i = 0; i < exponent; i++)
            res = res * base;
        return res;
    }

    static inline float flog2(float x) {
        if (x <= 0)
            return -128;
        MAY_ALIAS(int) *px = (int *) (&x);        // store address of float as pointer to long int
        float e = (float) (((*px & 0x7F800000) >> 23) -
                           0x7f); // shift right by 23 bits and subtract 127 = 0x7f => exponent
        *px = ((*px & 0x007FFFFF) | 0x3f800000);  // set exponent to 127 (i.e., 0)
        x -= 1.0;                                 // and calculate x-1.0
        x *= (1.441740
              + x * (-0.7077702 +
                     x * (0.4123442 + x * (-0.1903190 + x * 0.0440047)))); // 5'th order polynomial approx. of log(1+x)
        return x + e;
    }
};

#endif // COMPAT_MATHUTIL_H
