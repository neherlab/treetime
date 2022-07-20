# CentOS with official yum install procedure

FROM centos:7

# Setup Intel MKL
RUN yum-config-manager --add-repo https://yum.repos.intel.com/mkl/setup/intel-mkl.repo \
 && rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB \
 && yum install -y intel-mkl-2020.0-088 gcc pkg-config openssl-devel \
 && rm -rf /var/cache/yum/* \
 && yum clean all

# Setup Rust
# From official setting in https://github.com/rust-lang/docker-rust
ARG RUST_VERSION
ENV RUSTUP_HOME=/usr/local/rustup
ENV CARGO_HOME=/usr/local/cargo
ENV PATH=/usr/local/cargo/bin:$PATH
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --no-modify-path --default-toolchain ${RUST_VERSION}
# this may concern security issue for production use, but this container is designed for development use.
RUN chmod -R a+w $RUSTUP_HOME $CARGO_HOME

WORKDIR /src
RUN chmod -R a+w /src
