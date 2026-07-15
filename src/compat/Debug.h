// Minimal drop-in replacement for mmseqs' Debug.h / EXIT, so the ported
// taxonomy code compiles without linking against the mmseqs framework.
#ifndef COMPAT_DEBUG_H
#define COMPAT_DEBUG_H

#include <iostream>
#include <cstdlib>

class Debug {
public:
    // Match the level ordering used by mmseqs (higher == more verbose).
    enum LogLevel { NOTHING = 0, ERROR = 1, WARNING = 2, INFO = 3 };

    explicit Debug(int level) : level(level) {}

    template<typename T>
    Debug& operator<<(const T& value) {
        stream() << value;
        return *this;
    }

    // Support stream manipulators such as std::endl / std::flush.
    Debug& operator<<(std::ostream& (*manip)(std::ostream&)) {
        stream() << manip;
        return *this;
    }

private:
    std::ostream& stream() { return level <= WARNING ? std::cerr : std::cout; }
    int level;
};

#ifndef EXIT
#define EXIT(exitCode) std::exit(exitCode)
#endif

#endif // COMPAT_DEBUG_H
