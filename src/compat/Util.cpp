#include "Util.h"

#include <cstring>
#include <cstdlib>

// Behaviour matches mmseqs' Util::split: consecutive separators are collapsed
// (strtok semantics), so empty fields are not produced.
std::vector<std::string> Util::split(const std::string &str, const std::string &sep) {
    std::vector<std::string> arr;

    char *cstr = strdup(str.c_str());
    const char *csep = sep.c_str();
    char *rest;
    char *current = strtok_r(cstr, csep, &rest);
    while (current != NULL) {
        arr.emplace_back(current);
        current = strtok_r(NULL, csep, &rest);
    }
    free(cstr);

    return arr;
}
