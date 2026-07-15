#include "FileUtil.h"

#include <sys/stat.h>

bool FileUtil::fileExists(const char *fileName) {
    struct stat st;
    return stat(fileName, &st) == 0;
}
