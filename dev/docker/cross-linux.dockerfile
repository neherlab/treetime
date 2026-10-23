# syntax=docker/dockerfile:1
# check=experimental=all
#
# Cross image of a Linux release target. clang and lld of LLVM 22 compile and
# link for the target and for the build machine, each on the sysroot, crt files,
# libgcc, and libstdc++ of a crosstool-NG GCC 16 toolchain: glibc 2.17 for the
# gnu targets and the build machine, musl 1.2.6 for the musl targets. QEMU runs
# the binaries of the target.
FROM debian:12.15-slim@sha256:7ff6ba31527b9e336162fa061746c366c4b2032f4d2cff2078bea43f78efe3f8

SHELL ["bash", "-euxo", "pipefail", "-c"]

ARG CROSS_ARCH_DEBIAN
ARG CROSS_COMPILE
ARG CROSS_COMPILE_UPPER
ARG CROSS_COMPILE_SAFE
ARG CROSS_RUNNER

ENV HOST_TUPLE="x86_64-unknown-linux-gnu"
ENV HOST_TUPLE_UPPER="X86_64_UNKNOWN_LINUX_GNU"
ENV HOST_TUPLE_SAFE="x86_64_unknown_linux_gnu"
ENV CROSS_COMPILE="${CROSS_COMPILE}"
ENV CROSS_COMPILE_UPPER="${CROSS_COMPILE_UPPER}"
ENV CROSS_COMPILE_SAFE="${CROSS_COMPILE_SAFE}"
ENV CROSS_RUNNER="${CROSS_RUNNER}"

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& dpkg --add-architecture "${CROSS_ARCH_DEBIAN}" \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  bzip2 \
  ca-certificates \
  curl \
  file \
  git \
  "libc6:${CROSS_ARCH_DEBIAN}" \
  make \
  pigz \
  pixz \
  pkg-config \
  qemu-user \
  unzip \
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

# The target and the build machine get one GCC toolchain each, the same one when
# the target is x86_64-unknown-linux-gnu.
ENV LLVM_DIR="/opt/llvm"
ENV HOST_GCC_DIR="/opt/gcc-${HOST_TUPLE}"
ENV CROSS_GCC_DIR="/opt/gcc-${CROSS_COMPILE}"
COPY dev/docker/files/install-clang-target /
RUN set -euxo pipefail >/dev/null \
&& /install-image-tools "llvm=${LLVM_DIR}" "gcc-${HOST_TUPLE}=${HOST_GCC_DIR}" "gcc-${CROSS_COMPILE}=${CROSS_GCC_DIR}" \
&& "${LLVM_DIR}/bin/clang" --version \
&& for triplet in "${HOST_TUPLE}" "${CROSS_COMPILE}"; do \
  test -d "/opt/gcc-${triplet}/lib/gcc/${triplet}/16.2.0" \
  && test -d "/opt/gcc-${triplet}/${triplet}/sysroot" \
  && /install-clang-target "${LLVM_DIR}" "${triplet}" "/opt/gcc-${triplet}" \
  || exit 1; \
done \
&& rm /install-clang-target

# OpenBLAS, liblzma, and libzstd of the target: the prebuilt archives of
# .config/mise.image.toml, which cargo links through pkg-config with the [env]
# table that the install writes to /.cargo/config.toml.
COPY dev/docker/files/install-image-libraries /
RUN set -euxo pipefail >/dev/null \
&& /install-image-libraries "${CROSS_COMPILE}" \
&& rm /install-image-libraries

ENV PATH="${LLVM_DIR}/bin:${PATH}"
ENV CROSS_SYSROOT="${CROSS_GCC_DIR}/${CROSS_COMPILE}/sysroot"

ENV CC_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/${CROSS_COMPILE}-clang"
ENV CXX_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/${CROSS_COMPILE}-clang++"
ENV AR_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-ar"
ENV NM_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-nm"
ENV RANLIB_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-ranlib"
ENV STRIP_${CROSS_COMPILE_SAFE}="${LLVM_DIR}/bin/llvm-strip"
ENV LDD_${CROSS_COMPILE_SAFE}="${CROSS_GCC_DIR}/bin/${CROSS_COMPILE}-ldd"
ENV OBJDUMP_${CROSS_COMPILE_SAFE}="${CROSS_GCC_DIR}/bin/${CROSS_COMPILE}-objdump"
ENV STRINGS_${CROSS_COMPILE_SAFE}="${CROSS_GCC_DIR}/bin/${CROSS_COMPILE}-strings"
ENV CARGO_TARGET_${CROSS_COMPILE_UPPER}_LINKER="${LLVM_DIR}/bin/${CROSS_COMPILE}-clang"
ENV CARGO_TARGET_${CROSS_COMPILE_UPPER}_RUNNER="${CROSS_RUNNER}"

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
