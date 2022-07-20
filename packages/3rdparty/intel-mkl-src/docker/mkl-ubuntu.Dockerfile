# Ubuntu with official apt install procedure

FROM ubuntu:18.04

RUN apt update \
 && apt install -y \
      apt-utils \
      curl \
      gnupg \
      libssl-dev \
      pkg-config \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Install MKL
# https://software.intel.com/content/www/us/en/develop/articles/installing-intel-free-libs-and-python-apt-repo.html
RUN curl -sfL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB | apt-key add -
RUN curl -sfL https://apt.repos.intel.com/setup/intelproducts.list -o /etc/apt/sources.list.d/intelproducts.list

RUN apt update \
 && apt install -y intel-mkl-2020.0.088 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

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
