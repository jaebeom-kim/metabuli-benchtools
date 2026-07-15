cmake_minimum_required(VERSION 3.10)
project(metabuli-benchtools CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O2")

add_executable(benchtools
    src/main.cpp
    src/Parameters.cpp
    src/compat/Util.cpp
    src/compat/FileUtil.cpp
    src/taxonomy/NcbiTaxonomy.cpp
    src/taxonomy/TaxonomyWrapper.cpp
    src/tools/grade.cpp
    src/tools/makeBenchmarkSet.cpp
    src/tools/sampleQueries.cpp
)

target_include_directories(benchtools PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/compat
    ${CMAKE_SOURCE_DIR}/src/taxonomy
)

# OpenMP is optional: grade() parallelizes over input files when available.
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(benchtools PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(benchtools PRIVATE OPENMP=1)
endif()

install(TARGETS benchtools RUNTIME DESTINATION bin)
