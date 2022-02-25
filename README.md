# xAtlas

## Overview

xAtlas is a fast and scalable small variant caller that has been developed
at the Baylor College of Medicine Human Genome Sequencing Center.

## Building

xAtlas requires the following software components to build:

- [CMake](https://cmake.org/) (version 3.0 or higher required)
- [HTSlib](https://github.com/samtools/htslib) (version 1.6 or higher
recommended)
- libpthread (included with most Unix-like systems)

xAtlas can be built using CMake by running the following commands from the base
directory:

    $ mkdir build
    $ cd build
    $ cmake ..
    $ make

If HTSlib is installed in a location outside of CMake's library search path,
the HTSlib install prefix can be set for the CMake command by defining the
CMake variable `HTSLIB_PREFIX`:

    $ cmake -DHTSLIB_PREFIX=/path/to/htslib-1.x ..

## Usage

For example usage, see [GettingStarted](doc/GettingStarted.md).
