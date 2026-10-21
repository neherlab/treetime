# syntax=docker/dockerfile:1
# check=experimental=all
#
# Cross image of a macOS release target. The osxcross wrappers of the llvm
# flavor compile and link for the target with clang and ld64.lld of LLVM 22, the
# macOS 11.1 SDK, and compiler-rt, so the binaries load only system libraries.
# clang compiles and links for the build machine on the glibc 2.17 sysroot of a
# crosstool-NG GCC 16 toolchain.
FROM debian:12.15-slim@sha256:7ff6ba31527b9e336162fa061746c366c4b2032f4d2cff2078bea43f78efe3f8

SHELL ["bash", "-euxo", "pipefail", "-c"]

ARG CROSS_COMPILE
ARG CROSS_COMPILE_UPPER
ARG CROSS_COMPILE_SAFE
ARG CROSS_APPLE_TRIPLET
ARG CROSS_MACOS_DEPLOYMENT_TARGET

ENV HOST_TUPLE="x86_64-unknown-linux-gnu"
ENV HOST_TUPLE_UPPER="X86_64_UNKNOWN_LINUX_GNU"
ENV CROSS_COMPILE="${CROSS_COMPILE}"
ENV CROSS_COMPILE_UPPER="${CROSS_COMPILE_UPPER}"
ENV CROSS_COMPILE_SAFE="${CROSS_COMPILE_SAFE}"
ENV CROSS_APPLE_TRIPLET="${CROSS_APPLE_TRIPLET}"

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
  xz-utils \
  zstd \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

COPY dev/docker/files/fetch dev/docker/files/checksums /

ENV LLVM_DIR="/opt/llvm"
ENV HOST_GCC_DIR="/opt/gcc-${HOST_TUPLE}"
ENV OSX_CROSS_PATH="/opt/osxcross"
ENV PREFIX_CROSS="/opt/cross-${CROSS_COMPILE}"
COPY dev/docker/files/install-llvm dev/docker/files/install-gcc-cross dev/docker/files/install-clang-target dev/docker/files/install-osxcross dev/docker/files/install-openblas dev/docker/files/install-compression-libs /
RUN set -euxo pipefail >/dev/null \
&& /install-llvm "${LLVM_DIR}" \
&& /install-gcc-cross "${HOST_TUPLE}" "${HOST_GCC_DIR}" \
&& /install-clang-target "${LLVM_DIR}" "${HOST_TUPLE}" "${HOST_GCC_DIR}" \
&& /install-osxcross "${OSX_CROSS_PATH}" "${LLVM_DIR}" \
&& /install-openblas "${CROSS_COMPILE}" "/opt/openblas" \
&& /install-compression-libs "${CROSS_COMPILE}" "${PREFIX_CROSS}" \
&& rm /install-llvm /install-gcc-cross /install-clang-target /install-osxcross /install-openblas /install-compression-libs

ENV MACOSX_DEPLOYMENT_TARGET="${CROSS_MACOS_DEPLOYMENT_TARGET}"
ENV PATH="${OSX_CROSS_PATH}/bin:${LLVM_DIR}/bin:${PATH}"
ENV CROSS_SYSROOT="${OSX_CROSS_PATH}/SDK/MacOSX11.1.sdk"

ENV CC_${CROSS_COMPILE_SAFE}="${OSX_CROSS_PATH}/bin/${CROSS_APPLE_TRIPLET}-clang"
ENV CXX_${CROSS_COMPILE_SAFE}="${OSX_CROSS_PATH}/bin/${CROSS_APPLE_TRIPLET}-clang++"
ENV AR_${CROSS_COMPILE_SAFE}="${OSX_CROSS_PATH}/bin/${CROSS_APPLE_TRIPLET}-ar"
ENV NM_${CROSS_COMPILE_SAFE}="${OSX_CROSS_PATH}/bin/${CROSS_APPLE_TRIPLET}-nm"
ENV RANLIB_${CROSS_COMPILE_SAFE}="${OSX_CROSS_PATH}/bin/${CROSS_APPLE_TRIPLET}-ranlib"
ENV STRIP_${CROSS_COMPILE_SAFE}="${OSX_CROSS_PATH}/bin/${CROSS_APPLE_TRIPLET}-strip"
ENV OTOOL_${CROSS_COMPILE_SAFE}="${OSX_CROSS_PATH}/bin/${CROSS_APPLE_TRIPLET}-otool"
ENV CARGO_TARGET_${CROSS_COMPILE_UPPER}_LINKER="${OSX_CROSS_PATH}/bin/${CROSS_APPLE_TRIPLET}-clang"

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
ENV MISE_DATA_DIR="/opt/mise"
ENV MISE_CACHE_DIR="/tmp/mise/cache"
ENV MISE_STATE_DIR="/tmp/mise/state"
COPY dev/docker/files/install-mise dev/docker/files/mise-version /
COPY .config/mise.toml .config/mise.lock /tmp/mise/project/
RUN set -euxo pipefail >/dev/null \
&& /install-mise "/usr/local/bin" \
&& export MISE_TRUSTED_CONFIG_PATHS="/tmp/mise/project" \
&& mise -C "/tmp/mise/project" install "github:kunobi-ninja/kache" "protoc" \
&& ln -sf -t "/usr/local/bin" \
  "$(mise -C "/tmp/mise/project" which kache)" \
  "$(mise -C "/tmp/mise/project" which protoc)" \
&& chmod -R a+rX "${MISE_DATA_DIR}" \
&& rm -rf "/tmp/mise" "/install-mise" "/mise-version" "/usr/local/bin/mise" \
&& kache --version \
&& protoc --version
