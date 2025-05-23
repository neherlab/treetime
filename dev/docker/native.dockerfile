# syntax=docker/dockerfile:1
# check=experimental=all
FROM debian:12.8

SHELL ["bash", "-euxo", "pipefail", "-c"]

ENV HOST_TUPLE_DEBIAN="x86_64-linux-gnu"
ENV HOST_TUPLE="x86_64-unknown-linux-gnu"

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  bash \
  bash-completion \
  ca-certificates \
  curl \
  file \
  git \
  libc6-dev \
  lsb-release \
  make \
  parallel \
  pigz \
  pixz \
  pkg-config \
  python3 \
  python3-pip \
  sudo \
  tar \
  time \
  unzip \
  util-linux \
  xz-utils \
  zstd \
>/dev/null \
&& rm -rf /var/lib/apt/lists/* \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null



ENV HOST_GCC_DIR="/usr/local"
ENV HOSTCC="${HOST_GCC_DIR}/bin/gcc"
ENV HOSTCXX="${HOST_GCC_DIR}/bin/g++"
ENV HOSTFC="${HOST_GCC_DIR}/bin/gfortran"
ENV C_INCLUDE_PATH="/usr/include::/usr/local/include:/usr/include/${HOST_TUPLE_DEBIAN}"
ENV CPLUS_INCLUDE_PATH="${C_INCLUDE_PATH}"
ENV LIBRARY_PATH="/usr/lib:/usr/lib64:/usr/local/lib:/usr/local/lib64:/usr/lib/${HOST_TUPLE_DEBIAN}"
ENV LD_LIBRARY_PATH="/usr/lib:/usr/lib64:/usr/local/lib:/usr/local/lib64:/usr/lib/${HOST_TUPLE_DEBIAN}"

COPY --link "dev/docker/files/install-gcc" "/"
RUN /install-gcc "${HOST_GCC_DIR}"

COPY --link "dev/docker/files/install-llvm" "/"
RUN /install-llvm

COPY --link "dev/docker/files/install-protobuf" "/"
RUN /install-protobuf

COPY --link "dev/docker/files/install-hyperfine" "/"
RUN /install-hyperfine


ENV HOST_PREFIX="/usr"
ENV PKG_CONFIG_PATH="/usr/local/lib/pkgconfig:${HOST_PREFIX}/lib/pkgconfig"
ENV OPENBLAS_LIB_DIR="${HOST_PREFIX}/lib"
COPY --link "dev/docker/files/install-openblas" "/"
RUN /install-openblas "${HOST_TUPLE}" "${HOST_PREFIX}"

COPY --link "dev/docker/files/install-libbzip2" "/"
RUN /install-libbzip2 "${HOST_TUPLE}" "${HOST_PREFIX}"

COPY --link "dev/docker/files/install-liblzma" "/"
RUN /install-liblzma "${HOST_TUPLE}" "${HOST_PREFIX}"

COPY --link "dev/docker/files/install-libz" "/"
RUN /install-libz "${HOST_TUPLE}" "${HOST_PREFIX}"

COPY --link "dev/docker/files/install-libzstd" "/"
RUN /install-libzstd "${HOST_TUPLE}" "${HOST_PREFIX}"
ENV ZSTD_SYS_USE_PKG_CONFIG="1"
ENV LIBZ_SYS_STATIC="1"



ARG USER=user
ARG GROUP=user
ARG UID
ARG GID

ENV USER=$USER
ENV GROUP=$GROUP
ENV UID=$UID
ENV GID=$GID
ENV TERM="xterm-256color"
ENV HOME="/home/${USER}"

COPY --link "dev/docker/files/create-user" "/"
RUN /create-user


USER ${USER}


ENV CARGO_HOME="${HOME}/.cargo"
ENV PATH="${CARGO_HOME}/bin:${PATH}"
COPY --link --chown="${UID}:${GID}" "rust-toolchain.toml" "${CARGO_HOME}/rust-toolchain.toml"
COPY --link "dev/docker/files/install-rust" "/"
RUN set -euxo pipefail >/dev/null \
&& /install-rust "${HOST_TUPLE}" "${CARGO_HOME}"
