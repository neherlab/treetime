# syntax=docker/dockerfile:1
# check=experimental=all
#
# Cross image of a Windows release target (*-pc-windows-gnullvm). The
# llvm-mingw toolchain of LLVM 23 compiles and links for the target: clang and
# lld on the static mingw-w64 sysroot with the Universal C Runtime (UCRT),
# compiler-rt, and libunwind, so the binaries import only Windows system DLLs.
# clang compiles and links for the build machine on the glibc 2.17 sysroot of a
# crosstool-NG GCC 16 toolchain. Wine runs the x86_64 binaries.
FROM debian:13.7-slim@sha256:62bd423aeee912639ca8ee681c2b762cb3e25137e0a92b20be455898a9f034ad

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

# The <triple>-clang drivers of llvm-mingw load the sysroot, compiler-rt,
# libunwind, and lld from the LLVM prefix, so llvm-mingw is copied over the
# LLVM install.
ENV LLVM_DIR="/opt/llvm"
ENV HOST_GCC_DIR="/opt/gcc-${HOST_TUPLE}"
COPY dev/docker/files/install-clang-target /
RUN set -euxo pipefail >/dev/null \
&& /install-image-tools "llvm=${LLVM_DIR}" "llvm-mingw=/opt/llvm-mingw" "gcc-${HOST_TUPLE}=${HOST_GCC_DIR}" \
&& cp -a "/opt/llvm-mingw/." "${LLVM_DIR}/" \
&& "${LLVM_DIR}/bin/clang" --version \
&& for triple in x86_64-w64-mingw32 aarch64-w64-mingw32; do \
  "${LLVM_DIR}/bin/${triple}-clang" --version || exit 1; \
done \
&& test -d "${HOST_GCC_DIR}/lib/gcc/${HOST_TUPLE}/16.2.0" \
&& test -d "${HOST_GCC_DIR}/${HOST_TUPLE}/sysroot" \
&& /install-clang-target "${LLVM_DIR}" "${HOST_TUPLE}" "${HOST_GCC_DIR}" \
&& rm /install-clang-target

# OpenBLAS, liblzma, and libzstd of the target: the prebuilt archives of
# .config/mise.image.toml, which cargo links through pkg-config with the [env]
# table that the install writes to /.cargo/config.toml. The install fails when
# an archive was built with another toolchain or CPU flags than this image uses.
COPY dev/docker/files/install-image-libraries dev/lib/dist-flags.sh /
RUN set -euxo pipefail >/dev/null \
&& /install-image-libraries "${CROSS_COMPILE}" \
&& rm /install-image-libraries /dist-flags.sh

ENV PATH="${LLVM_DIR}/bin:${PATH}"

ENV CC_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/${CROSS_MINGW_TRIPLET}-clang"
ENV CXX_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/${CROSS_MINGW_TRIPLET}-clang++"
ENV AR_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-ar"
ENV NM_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-nm"
ENV RANLIB_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-ranlib"
ENV STRIP_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-strip"
ENV OBJDUMP_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-objdump"
ENV CARGO_TARGET_${CROSS_COMPILE_UPPER}_LINKER="${LLVM_DIR}/bin/${CROSS_MINGW_TRIPLET}-clang"

ENV HOST_CC="${LLVM_DIR}/bin/${HOST_TUPLE}-clang"
ENV HOST_CXX="${LLVM_DIR}/bin/${HOST_TUPLE}-clang++"
ENV HOST_AR="${LLVM_DIR}/bin/llvm-ar"
ENV CARGO_TARGET_${HOST_TUPLE_UPPER}_LINKER="${LLVM_DIR}/bin/${HOST_TUPLE}-clang"

# The Rust toolchain of rust-toolchain.toml with the target. dev/docker/run runs
# as the host user with HOME=/tmp/home and the cargo home in the checkout.
ENV RUSTUP_HOME="/usr/local/rustup"
ENV CARGO_HOME="/usr/local/cargo"
ENV PATH="/usr/local/cargo/bin:${PATH}"
COPY dev/docker/files/install-rust /
COPY rust-toolchain.toml /tmp/rust/
RUN set -euxo pipefail >/dev/null \
&& /install-rust "/tmp/rust" "${CROSS_COMPILE}" \
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
