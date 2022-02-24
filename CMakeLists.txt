cmake_minimum_required(VERSION 3.0 FATAL_ERROR)

# target: xatlas

project(xatlas
    VERSION 0.3
    LANGUAGES CXX
)

add_compile_options(-Wall -Wextra -pedantic -O2)
#add_compile_options(-Wall -Wextra -pedantic -O0 -g)

# xatlas: sources and executable

include_directories(include)

set(XATLAS_SOURCES
    src/Bam.cc
    src/CoverageCounter.cc
    src/EventScanner.cc
    src/GvcfBlock.cc
    src/IndelEvent.cc
    src/Logit.cc
    src/ReferenceSequence.cc
    src/SnpEvent.cc
    src/VcfWriter.cc
    src/Xatlas.cc
)
add_executable(xatlas ${XATLAS_SOURCES})

set_target_properties(xatlas
    PROPERTIES
    CXX_STANDARD 11
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
)

# xatlas: dependencies

## HTSlib

#set(HTSLIB_PREFIX "/hgsc_software/htslib/htslib-1.10")
if(DEFINED HTSLIB_PREFIX)
    # using custom htslib prefix
    find_library(HTSLIB
        NAMES hts
        PATHS "${HTSLIB_PREFIX}/lib"
        REQUIRED
    )
    target_include_directories(xatlas
        PUBLIC
        "${HTSLIB_PREFIX}/include"
    )
else()
    # using default search paths
    find_library(HTSLIB
        NAMES hts
        REQUIRED
    )
endif()

## pthread

if(DEFINED PTHREADLIB_PREFIX)
    # using custom htslib prefix
    find_library(PTHREADLIB
        NAMES pthread
        PATHS "${PTHREADLIB_PREFIX}/lib"
        REQUIRED
    )
    target_include_directories(xatlas
        PUBLIC
        "${PTHREADLIB_PREFIX}/include")
else()
    # using default search paths
    find_library(PTHREADLIB
        NAMES pthread
        REQUIRED
    )
endif()

## set linked libraries

set(LIBRARIES ${PTHREADLIB} ${HTSLIB})
target_link_libraries(xatlas PUBLIC ${LIBRARIES})

# xatlas: install targets

install(TARGETS xatlas DESTINATION bin)
