FROM ubuntu:20.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq update && apt-get -yq install cmake build-essential wget libcurl4-openssl-dev libz-dev liblzma-dev libbz2-dev && \
    wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 -O htslib.tar.bz2 && \
    tar -xjvf htslib.tar.bz2 && \
    cd htslib-1.15.1 && \
    make && \
    make install 

ADD . /opt/xatlas-source
WORKDIR /opt/xatlas-source
RUN mkdir build &&  \
    cd build  &&  \
    cmake .. && \
    make && \
    mv xatlas /bin/

WORKDIR /data
ENTRYPOINT ["xatlas"]
