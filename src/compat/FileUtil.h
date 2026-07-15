// Minimal subset of mmseqs' FileUtil.h. The ported taxonomy code only ever
// calls FileUtil::fileExists(const char*).
#ifndef COMPAT_FILEUTIL_H
#define COMPAT_FILEUTIL_H

class FileUtil {
public:
    static bool fileExists(const char *fileName);
};

#endif // COMPAT_FILEUTIL_H
