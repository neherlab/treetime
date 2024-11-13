# syntax=docker/dockerfile:1
# check=experimental=all
FROM debian:12.7

SHELL ["bash", "-euxo", "pipefail", "-c"]

ENV HOST_TUPLE_DEBIAN="x86_64-linux-gnu"
ENV HOST_TUPLE="x86_64-unknown-linux-gnu"

ARG CROSS_ARCH_DEBIAN
ARG CROSS_ARCH
ARG CROSS_COMPILE
ARG CROSS_COMPILE_UPPER
ARG CROSS_GCC_TRIPLET
ARG CROSS_RUNNER

ENV HOST_TUPLE="x86_64-linux-gnu"
ENV CROSS_ARCH_DEBIAN="${CROSS_ARCH_DEBIAN}"
ENV CROSS_ARCH="${CROSS_ARCH}"
ENV CROSS_COMPILE="${CROSS_COMPILE}"
ENV CROSS_COMPILE_UPPER="${CROSS_COMPILE_UPPER}"
ENV CROSS_GCC_TRIPLET="${CROSS_GCC_TRIPLET}"
ENV CROSS_RUNNER="${CROSS_RUNNER}"

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& dpkg --add-architecture ${CROSS_ARCH_DEBIAN} \
&& apt-get update -qq --yes \
&& apt-get install -qq --no-install-recommends --yes \
  bash \
  ca-certificates \
  curl \
  file \
  git \
  libc6-dev \
  make \
  pigz \
  pixz \
  pkg-config \
  sudo \
  tar \
  time \
  unzip \
  util-linux \
  xz-utils \
  zstd \
>/dev/null \
&& if [[ "${CROSS_COMPILE}" =~ (linux) ]]; then apt-get install -qq --no-install-recommends --yes  \
  libc6:${CROSS_ARCH_DEBIAN} \
  qemu-user \
>/dev/null \
;fi \
&& if [[ "${CROSS_COMPILE}" =~ (mingw|windows) ]]; then apt-get install -qq --no-install-recommends --yes  \
  wine64 \
>/dev/null \
;fi \
&& rm -rf /var/lib/apt/lists/* \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null



ENV HOST_GCC_DIR="/usr/local"
ENV HOSTCC="${HOST_GCC_DIR}/bin/gcc"
ENV HOSTCXX="${HOST_GCC_DIR}/bin/g++"
ENV HOSTFC="${HOST_GCC_DIR}/bin/gfortran"
ENV LIBRARY_PATH="/usr/lib:/usr/lib64:/usr/local/lib:/usr/local/lib64:/usr/lib/${HOST_TUPLE_DEBIAN}"
ENV LD_LIBRARY_PATH="/usr/lib:/usr/lib64:/usr/local/lib:/usr/local/lib64:/usr/lib/${HOST_TUPLE_DEBIAN}"

COPY --link "dev/docker/files/install-gcc" "/"
RUN /install-gcc "${HOST_GCC_DIR}"

COPY --link "dev/docker/files/install-llvm" "/"
RUN /install-llvm

COPY --link "dev/docker/files/install-protobuf" "/"
RUN /install-protobuf


ENV CROSS_GCC_DIR="/opt/gcc-${CROSS_GCC_TRIPLET}"
ENV PATH="${CROSS_GCC_DIR}/bin:${PATH}"
ENV CROSS_SYSROOT="${CROSS_GCC_DIR}/${CROSS_GCC_TRIPLET}/sysroot"
ENV CC_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-gcc"
ENV CXX_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-g++"
ENV FC_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-gfortran"
ENV ADDR2LINE_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-addr2line"
ENV AR_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-gcc-ar"
ENV AS_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-as"
ENV CPP_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-cpp"
ENV ELFEDIT_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-elfedit"
ENV LD_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-ld"
ENV LDD_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-ldd"
ENV NM_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-gcc-nm"
ENV OBJCOPY_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-objcopy"
ENV OBJDUMP_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-objdump"
ENV RANLIB_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-gcc-ranlib"
ENV READELF_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-readelf"
ENV SIZE_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-size"
ENV STRINGS_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-strings"
ENV STRIP_${CROSS_COMPILE}="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-strip"
ENV BINDGEN_EXTRA_CLANG_ARGS_${CROSS_COMPILE}="--sysroot=${CROSS_SYSROOT}"
ENV CARGO_TARGET_${CROSS_COMPILE_UPPER}_AR="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-gcc-ar"
ENV CARGO_TARGET_${CROSS_COMPILE_UPPER}_LINKER="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-gcc"
ENV CARGO_TARGET_${CROSS_COMPILE_UPPER}_STRIP="${CROSS_GCC_DIR}/bin/${CROSS_GCC_TRIPLET}-strip"
ENV CARGO_TARGET_${CROSS_COMPILE_UPPER}_RUNNER="${CROSS_RUNNER}"

COPY --link "dev/docker/files/install-gcc-cross" "/"
RUN /install-gcc-cross "${CROSS_GCC_TRIPLET}" "${CROSS_GCC_DIR}"


ENV PREFIX_CROSS="/usr/local/${CROSS_COMPILE}"
ENV PKG_CONFIG_ALLOW_CROSS="1"
ENV PKG_CONFIG_SYSROOT_DIR="${PREFIX_CROSS}"
ENV PKG_CONFIG_PATH="${PREFIX_CROSS}/lib/pkgconfig:/usr/local/lib/pkgconfig"


ENV OPENBLAS_LIB_DIR="${PREFIX_CROSS}/lib"
COPY --link "dev/docker/files/install-openblas" "/"
RUN /install-openblas "${CROSS_COMPILE}" "${PREFIX_CROSS}"

COPY --link "dev/docker/files/install-libbzip2" "/"
RUN /install-libbzip2 "${CROSS_COMPILE}" "${PREFIX_CROSS}"

COPY --link "dev/docker/files/install-liblzma" "/"
RUN /install-liblzma "${CROSS_COMPILE}" "${PREFIX_CROSS}"

COPY --link "dev/docker/files/install-libz" "/"
RUN /install-libz "${CROSS_COMPILE}" "${PREFIX_CROSS}"

COPY --link "dev/docker/files/install-libzstd" "/"
RUN /install-libzstd "${CROSS_COMPILE}" "${PREFIX_CROSS}"
ENV ZSTD_SYS_USE_PKG_CONFIG="1"


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
&& /install-rust "${CROSS_COMPILE}" "${CARGO_HOME}"



ENV WINEPREFIX="${HOME}/wine"
ENV WINEARCH="win64"
ENV WINEDEBUG="err+all,err-winediag,err-ole,fixme-all"
COPY --link "dev/docker/files/install-bryptprimitives" "/"
RUN if [[ "${CROSS_COMPILE}" =~ (mingw|windows) ]]; then /install-bryptprimitives "${WINEPREFIX}/drive_c/windows/system32"; fi
