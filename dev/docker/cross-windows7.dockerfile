# syntax=docker/dockerfile:1
# check=experimental=all
#
# Cross image of the Windows 7 target x86_64-win7-windows-gnu. Rust has no
# prebuilt standard library for this tier 3 target, so cargo builds it from the
# rust-src component (-Zbuild-std, dev/cross/build). The crosstool-NG GCC 16
# MinGW toolchain on msvcrt compiles and links for the target, because Windows 7
# has no Universal C Runtime without an update; OpenBLAS is the msvcrt build of
# binarylandia/build_openblas with the static Fortran runtime of flang. clang
# compiles and links for the build machine on the glibc 2.17 sysroot of a
# crosstool-NG GCC 16 toolchain. Wine runs the binaries.
FROM debian:12.15-slim@sha256:7ff6ba31527b9e336162fa061746c366c4b2032f4d2cff2078bea43f78efe3f8

SHELL ["bash", "-euxo", "pipefail", "-c"]

ARG CROSS_COMPILE
ARG CROSS_COMPILE_UPPER
ARG CROSS_COMPILE_SAFE
ARG CROSS_MINGW_TRIPLET
ARG CROSS_RUNNER

ENV HOST_TUPLE="x86_64-unknown-linux-gnu"
ENV HOST_TUPLE_UPPER="X86_64_UNKNOWN_LINUX_GNU"
ENV CROSS_COMPILE="${CROSS_COMPILE}"
ENV CROSS_COMPILE_UPPER="${CROSS_COMPILE_UPPER}"
ENV CROSS_COMPILE_SAFE="${CROSS_COMPILE_SAFE}"
ENV CROSS_MINGW_TRIPLET="${CROSS_MINGW_TRIPLET}"
ENV CROSS_RUNNER="${CROSS_RUNNER}"

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  ca-certificates \
  curl \
  file \
  git \
  make \
  parallel \
  pigz \
  pixz \
  pkg-config \
  unzip \
  wine64 \
  xz-utils \
  zstd \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

COPY dev/docker/files/fetch dev/docker/files/checksums /

# mise installs the build toolchains of .config/mise.image.toml here and the
# tools of .config/mise.toml further down, each at the URLs and checksums of its
# lockfile, so a change of the development tools keeps the toolchain layers.
ENV MISE_DATA_DIR="/opt/mise"
ENV MISE_CACHE_DIR="/tmp/mise/cache"
ENV MISE_STATE_DIR="/tmp/mise/state"
COPY dev/docker/files/install-mise dev/docker/files/mise-version /
RUN set -euxo pipefail >/dev/null \
&& /install-mise "/usr/local/bin" \
&& rm /install-mise /mise-version
COPY dev/docker/files/install-image-tools /
COPY .config/mise.image.toml /mise-image/mise.toml
COPY .config/mise.image.lock /mise-image/mise.lock

ENV LLVM_DIR="/opt/llvm"
ENV HOST_GCC_DIR="/opt/gcc-${HOST_TUPLE}"
ENV CROSS_GCC_DIR="/opt/gcc-${CROSS_MINGW_TRIPLET}"
COPY dev/docker/files/install-clang-target dev/docker/files/install-openblas /
RUN set -euxo pipefail >/dev/null \
&& /install-image-tools "llvm=${LLVM_DIR}" "gcc-${HOST_TUPLE}=${HOST_GCC_DIR}" "gcc-${CROSS_MINGW_TRIPLET}=${CROSS_GCC_DIR}" \
&& "${LLVM_DIR}/bin/clang" --version \
&& for triplet in "${HOST_TUPLE}" "${CROSS_MINGW_TRIPLET}"; do \
  test -d "/opt/gcc-${triplet}/lib/gcc/${triplet}/16.2.0" \
  && test -d "/opt/gcc-${triplet}/${triplet}/sysroot" \
  || exit 1; \
done \
&& /install-clang-target "${LLVM_DIR}" "${HOST_TUPLE}" "${HOST_GCC_DIR}" \
&& /install-openblas "${CROSS_COMPILE}" "/opt/openblas" \
&& rm /install-clang-target /install-openblas

ENV PATH="${LLVM_DIR}/bin:${CROSS_GCC_DIR}/bin:${PATH}"

ENV CC_${CROSS_COMPILE_SAFE}="${CROSS_GCC_DIR}/bin/${CROSS_MINGW_TRIPLET}-gcc"
ENV CXX_${CROSS_COMPILE_SAFE}="${CROSS_GCC_DIR}/bin/${CROSS_MINGW_TRIPLET}-g++"
ENV AR_${CROSS_COMPILE_SAFE}="${CROSS_GCC_DIR}/bin/${CROSS_MINGW_TRIPLET}-ar"
ENV RANLIB_${CROSS_COMPILE_SAFE}="${CROSS_GCC_DIR}/bin/${CROSS_MINGW_TRIPLET}-ranlib"
ENV OBJDUMP_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-objdump"
ENV READOBJ_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-readobj"
ENV CARGO_TARGET_${CROSS_COMPILE_UPPER}_LINKER="${CROSS_GCC_DIR}/bin/${CROSS_MINGW_TRIPLET}-gcc"

ENV HOST_CC="${LLVM_DIR}/bin/${HOST_TUPLE}-clang"
ENV HOST_CXX="${LLVM_DIR}/bin/${HOST_TUPLE}-clang++"
ENV HOST_AR="${LLVM_DIR}/bin/llvm-ar"
ENV CARGO_TARGET_${HOST_TUPLE_UPPER}_LINKER="${LLVM_DIR}/bin/${HOST_TUPLE}-clang"

# The functions that Windows 7 lacks, from the mingw-w64 headers of the
# toolchain, for the Windows 7 import check of dev/cross/check.
ENV WINDOWS8_FUNCTIONS="/opt/windows7/windows8-functions"
COPY dev/docker/files/write-windows8-functions /
RUN set -euxo pipefail >/dev/null \
&& mkdir -p "$(dirname "${WINDOWS8_FUNCTIONS}")" \
&& CC="${LLVM_DIR}/bin/clang --target=${CROSS_MINGW_TRIPLET} --sysroot=${CROSS_GCC_DIR}/${CROSS_MINGW_TRIPLET}/sysroot" \
  /write-windows8-functions "${WINDOWS8_FUNCTIONS}" \
&& rm /write-windows8-functions

# The Rust toolchain of rust-toolchain.toml with the source of the standard
# library, which cargo compiles for this target. dev/docker/run runs as the host
# user with HOME=/tmp/home and the cargo home in the checkout.
ENV RUSTUP_HOME="/usr/local/rustup"
ENV CARGO_HOME="/usr/local/cargo"
ENV PATH="/usr/local/cargo/bin:${PATH}"
COPY dev/docker/files/install-rust /
COPY rust-toolchain.toml /tmp/rust/
RUN set -euxo pipefail >/dev/null \
&& /install-rust "/tmp/rust" \
&& rustup component add rust-src \
&& chmod -R a+w "${RUSTUP_HOME}" "${CARGO_HOME}" \
&& rm -rf /install-rust /tmp/rust \
&& mkdir -p "/tmp/home" \
&& chmod 1777 "/tmp/home"

# Windows binaries run under Wine; dev/cross/run creates the Wine prefix.
ENV WINEARCH="win64"
ENV WINEDEBUG="err+all,err-winediag,err-ole,fixme-all"

# kache, the optional compiler cache of dev/cross/build, and protoc, which
# prost-build runs for the UShER protobuf schemas, at the versions, URLs, and
# sha256 of .config/mise.lock, so the cross and development images run the same
# tools.
COPY .config/mise.toml .config/mise.lock /tmp/mise/project/
RUN set -euxo pipefail >/dev/null \
&& export MISE_TRUSTED_CONFIG_PATHS="/tmp/mise/project" \
&& mise -C "/tmp/mise/project" install "github:kunobi-ninja/kache" "protoc" \
&& ln -sf -t "/usr/local/bin" \
  "$(mise -C "/tmp/mise/project" which kache)" \
  "$(mise -C "/tmp/mise/project" which protoc)" \
&& chmod -R a+rX "${MISE_DATA_DIR}" \
&& rm -rf "/tmp/mise" "/install-image-tools" "/mise-image" "/usr/local/bin/mise" \
&& kache --version \
&& protoc --version
