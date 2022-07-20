# rustmath/mkl-rust

ARG RUST_VERSION
FROM rust:${RUST_VERSION}

# Setup Intel-MKL
WORKDIR /mkl

RUN apt update \
 && apt install -y cpio \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

COPY ./silent.cfg /mkl/
RUN curl -sfLO http://registrationcenter-download.intel.com/akdlm/irc_nas/tec/16533/l_mkl_2020.1.217.tgz \
 && tar xf l_mkl_2020.1.217.tgz   \
 && cd l_mkl_2020.1.217           \
 && ./install.sh -s ../silent.cfg \
 && rm -rf /mkl

# Setup linker to find shared library
COPY intel-mkl.conf /etc/ld.so.conf.d/
RUN ldconfig

# Setup pkg-config
ENV PKG_CONFIG_PATH /opt/intel/mkl/bin/pkgconfig
RUN sed -i "s/MKLROOT/prefix/g" ${PKG_CONFIG_PATH}/*.pc

# Setup basic rust development tools
WORKDIR /src
RUN cargo install cargo-tarpaulin
