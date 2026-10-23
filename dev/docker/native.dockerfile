# syntax=docker/dockerfile:1
# check=experimental=all
#
# Development image: Rust and Bun toolchains, lint and test tools, and the
# libraries the desktop app needs. dev/docker/run builds it and runs every
# command as the host user, with HOME=/tmp/home and the cargo home in the
# checkout. C, C++, and Rust code for x86_64-unknown-linux-gnu compiles and
# links with the compiler, sysroot, and libraries of the shipped Linux build, so
# the binaries built here need glibc 2.17 of the GCC 16.2.0 sysroot, not the
# glibc of Ubuntu, and profiling builds also run on the host.
FROM ubuntu:noble-20260930.1@sha256:68ee22004a9109b44973313b02ff7b1de250f1de0f5ecb628d73f579aaff33ac
SHELL ["bash", "-euxo", "pipefail", "-c"]

# libatomic1 is a run-time library of mold.
RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  binutils \
  ca-certificates \
  curl \
  file \
  git \
  libatomic1 \
  libssl-dev \
  make \
  pigz \
  pixz \
  pkg-config \
  python3 \
  unzip \
  xz-utils \
  zstd \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

# Runtime libraries of Electron, for the desktop app in development.
RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  dbus \
  libasound2t64 \
  libatk-bridge2.0-0t64 \
  libatk1.0-0t64 \
  libatspi2.0-0t64 \
  libcairo2 \
  libcups2t64 \
  libdbus-1-3 \
  libdrm2 \
  libgbm1 \
  libgl1 \
  libgl1-mesa-dri \
  libglib2.0-0t64 \
  libgtk-3-0t64 \
  libnss3 \
  libpango-1.0-0 \
  libx11-6 \
  libxcb1 \
  libxcomposite1 \
  libxdamage1 \
  libxext6 \
  libxfixes3 \
  libxkbcommon0 \
  libxrandr2 \
  libxshmfence1 \
  xdg-utils \
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

# clang of LLVM 22 compiles and links x86_64-unknown-linux-gnu code on the
# sysroot (glibc 2.17), crt files, libgcc, and libstdc++ of the GCC 16.2.0
# toolchain, as in dev/docker/cross-linux.dockerfile. install-clang-target
# writes the configuration x86_64-unknown-linux-gnu.cfg next to clang, which
# clang reads for every call that targets this triple, its default, so `cc` and
# `c++` use it too, as do the cc crate, the build scripts, and rustc, whose
# default linker is `cc`. CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER stays
# unset: it would replace the dylint-link linker of the lint libraries
# (dev/lints/*/.cargo/config.toml). Development builds select mold in their
# rustflags; clang takes the last -fuse-ld, so the rustflags win over
# -fuse-ld=lld of the configuration.
ENV PATH="/opt/llvm/bin:${PATH}"
COPY dev/docker/files/install-clang-target /
RUN set -euxo pipefail >/dev/null \
&& /install-image-tools "llvm=/opt/llvm" "gcc-x86_64-unknown-linux-gnu=/opt/gcc-x86_64-unknown-linux-gnu" \
&& clang --version \
&& test -d "/opt/gcc-x86_64-unknown-linux-gnu/lib/gcc/x86_64-unknown-linux-gnu/16.2.0" \
&& test -d "/opt/gcc-x86_64-unknown-linux-gnu/x86_64-unknown-linux-gnu/sysroot" \
&& /install-clang-target "/opt/llvm" "x86_64-unknown-linux-gnu" "/opt/gcc-x86_64-unknown-linux-gnu" \
&& ln -sf "/opt/llvm/bin/x86_64-unknown-linux-gnu-clang" "/usr/local/bin/cc" \
&& ln -sf "/opt/llvm/bin/x86_64-unknown-linux-gnu-clang++" "/usr/local/bin/c++" \
&& printf 'int main(void) { return 0; }\n' | cc -x c -o /tmp/cc-check - \
&& printf '#include <string>\nint main() { return std::string("a").size() - 1; }\n' | c++ -x c++ -o /tmp/cxx-check - \
&& /tmp/cxx-check \
&& llvm-objdump -T /tmp/cc-check | grep -qE 'GLIBC_2\.2\.5\).*__libc_start_main' \
&& rm /install-clang-target /tmp/cc-check /tmp/cxx-check
ENV CC_x86_64_unknown_linux_gnu="/opt/llvm/bin/x86_64-unknown-linux-gnu-clang"
ENV CXX_x86_64_unknown_linux_gnu="/opt/llvm/bin/x86_64-unknown-linux-gnu-clang++"
ENV AR_x86_64_unknown_linux_gnu="/opt/llvm/bin/llvm-ar"
ENV HOST_CC="/opt/llvm/bin/x86_64-unknown-linux-gnu-clang"
ENV HOST_CXX="/opt/llvm/bin/x86_64-unknown-linux-gnu-clang++"
ENV HOST_AR="/opt/llvm/bin/llvm-ar"

# The Rust toolchain of rust-toolchain.toml, then the pinned nightly of the lint
# libraries, the toolchain that cargo-hawk is built against, and the nightly of
# the threaded WebAssembly build, which rebuilds std with atomics from rust-src
# (-Zbuild-std).
ENV RUSTUP_HOME="/usr/local/rustup"
ENV CARGO_HOME="/usr/local/cargo"
ENV PATH="/usr/local/cargo/bin:${PATH}"
COPY dev/docker/files/install-rust /
COPY rust-toolchain.toml /tmp/rust/
COPY dev/lints/dylint-custom/rust-toolchain.toml /tmp/lints/dylint-custom/
COPY dev/lints/dylint-mordant/rust-toolchain.toml /tmp/lints/dylint-mordant/
COPY dev/lints/dylint-trailofbits/rust-toolchain.toml /tmp/lints/dylint-trailofbits/
COPY dev/docker/files/hawk-toolchain dev/docker/files/wasm-threads-toolchain /tmp/
# hadolint ignore=DL3003
RUN set -euxo pipefail >/dev/null \
&& /install-rust "/tmp/rust" \
&& for dir in /tmp/lints/*; do (cd "${dir}" && rustup toolchain install) || exit 1; done \
&& rustup toolchain install "$(cat /tmp/hawk-toolchain)" --profile minimal --component rustc-dev,llvm-tools-preview,rust-src \
&& rustup toolchain install "$(cat /tmp/wasm-threads-toolchain)" --profile minimal --component rust-src,clippy \
&& chmod -R a+w "${RUSTUP_HOME}" "${CARGO_HOME}" \
&& rm -rf /install-rust /tmp/rust /tmp/lints /tmp/hawk-toolchain /tmp/wasm-threads-toolchain

# mise installs every tool of .config/mise.toml at the URL and sha256 of .config/mise.lock and
# links their executables into /usr/local/bin. The cargo registry that source
# builds fetch is a BuildKit cache, so it stays out of the image. The host
# libraries of .config/mise.toml stay disabled, in this layer and in every `just`
# recipe, because the image links its own prebuilt libraries (see below).
ENV MISE_DISABLE_TOOLS="openblas,liblzma"
COPY .config/mise.toml .config/mise.lock /tmp/mise/project/
RUN --mount=type=cache,target=/usr/local/cargo/registry,sharing=locked \
  set -euxo pipefail >/dev/null \
&& export MISE_TRUSTED_CONFIG_PATHS="/tmp/mise/project" \
&& mise -C "/tmp/mise/project" install \
&& bin_paths="$(mise -C "/tmp/mise/project" bin-paths)" \
&& printf '%s\n' "${bin_paths}" | while IFS= read -r dir; do \
  find -L "${dir}" -mindepth 1 -maxdepth 1 -type f -executable -exec ln -sf -t "/usr/local/bin" {} + || exit 1; \
done \
&& chmod -R a+rX "${MISE_DATA_DIR}" \
&& rm -rf "/tmp/mise" \
&& just --version \
&& cargo nextest --version \
&& cargo dylint --version

# Emscripten builds the WebAssembly back end of the web2 app for
# wasm32-unknown-emscripten: rustc links through emcc, which finds Node on PATH.
# Only the emcc front end goes on PATH, so the clang of Emscripten stays out of
# native builds. emcc runs the configuration
# .emscripten as Python at start: it names the LLVM and Binaryen of the archive
# and drops the include paths of the native GCC, which the clang of emcc would
# search too, mixing the glibc headers of the host into WebAssembly compiles.
# The cache holds the prebuilt system libraries, the threaded ones included,
# and stays writable for its lock and sanity files.
ENV PATH="/opt/emscripten/emscripten:${PATH}"
RUN set -euxo pipefail >/dev/null \
&& /install-image-tools "http:emscripten=/opt/emscripten" \
&& printf '%s\n' \
  "import os" \
  "" \
  "LLVM_ROOT = '/opt/emscripten/bin'" \
  "BINARYEN_ROOT = '/opt/emscripten'" \
  "" \
  "for name in ('C_INCLUDE_PATH', 'CPLUS_INCLUDE_PATH'):" \
  "  os.environ.pop(name, None)" \
  >"/opt/emscripten/emscripten/.emscripten" \
&& version="$(MISE_TRUSTED_CONFIG_PATHS="/mise-image" mise -C "/mise-image" current "http:emscripten")" \
&& grep -qxF -e "${version}-git" "/opt/emscripten/emscripten/emscripten-version.txt" \
&& emcc --version \
&& chmod -R a+rwX "/opt/emscripten/emscripten/cache"

# OpenBLAS, liblzma, and libzstd of the native and the WebAssembly builds: the
# prebuilt archives of .config/mise.image.toml, which cargo links through
# pkg-config with the [env] table that the install writes to
# /.cargo/config.toml, and for the threaded WebAssembly module with
# /.cargo/wasm-mt-env.toml. The install fails when an archive was built with
# another toolchain, CPU flags, or Emscripten version than this image uses.
COPY dev/docker/files/install-image-libraries dev/lib/dist-flags.sh /
RUN set -euxo pipefail >/dev/null \
&& /install-image-libraries "x86_64-unknown-linux-gnu" "wasm32-unknown-emscripten" \
&& rm -rf /install-image-libraries /dist-flags.sh /install-image-tools /mise-image "/tmp/mise"

# dev/docker/run runs as the host user with HOME=/tmp/home.
RUN set -euxo pipefail >/dev/null \
&& mkdir -p "/tmp/home" \
&& chmod 1777 "/tmp/home"
