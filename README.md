# xAtlas

## Overview

xAtlas is a fast and retrainable small variant caller that has been developed
at the Baylor College of Medicine Human Genome Sequencing Center.

## Building

xAtlas requires a recent version of HTSlib (versions 1.3 and higher are known
to work, version 1.4 or higher is suggested).

To build xAtlas, run:

    autoconf
    ./configure
    make

xAtlas multithreading support is enabled by running ``configure`` as:

    ./configure --enable-multithreading

## Usage

xAtlas supports two types of multithreading. The ``-P`` option splits the
processes for reading the input alignment, processing SNPs, and processing
indels into three separate threads. The ``-t`` option allows a variable number
of HTSlib decompression threads to be spawned.
