// Minimal subset of mmseqs' Util.h. Only Util::split is used by the ported
// code. EXIT lives in Debug.h; it is pulled in here because upstream code
// includes "Util.h" expecting EXIT to be available.
#ifndef COMPAT_UTIL_H
#define COMPAT_UTIL_H

#include <string>
#include <vector>
#include <sstream>
#include "Debug.h"

class Util {
public:
    static std::vector<std::string> split(const std::string &str, const std::string &sep);
};

// mmseqs exposes SSTR as a family of type-specialized stringify helpers. The
// ported taxonomy code only ever stringifies streamable scalars, so a single
// generic template is sufficient.
template<typename T>
inline std::string SSTR(const T &value) {
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

#endif // COMPAT_UTIL_H
