# syntax=docker/dockerfile:1
# check=experimental=all
#
# Development image: Rust and Bun toolchains, lint and test tools, and the
# libraries the desktop app needs. dev/docker/run builds it and runs every
# command as the host user, with HOME=/tmp/home and the cargo home in the
# checkout. Ubuntu 24.04 keeps the glibc of the binaries built here no newer than
# on current hosts, so profiling builds run on the host.
FROM ubuntu:noble-20260930.1@sha256:68ee22004a9109b44973313b02ff7b1de250f1de0f5ecb628d73f579aaff33ac
SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  binutils \
  ca-certificates \
  curl \
  file \
  git \
  libc6-dev \
  libfontconfig-dev \
  libgcc-13-dev \
  libssl-dev \
  libstdc++-13-dev \
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

# clang of LLVM 22 is the C compiler and link driver (`cc` and `c++`), on the C
# library and the GCC runtime of Ubuntu; the native builds link with mold.
# OpenBLAS, liblzma, and libzstd are the prebuilt archives of binarylandia, and
# /.cargo/config.toml, which their install scripts write, makes cargo link them
# instead of compiling them or asking for shared libraries.
ENV PATH="/opt/llvm/bin:${PATH}"
COPY dev/docker/files/install-openblas dev/docker/files/install-compression-libs /
RUN set -euxo pipefail >/dev/null \
&& /install-image-tools "llvm=/opt/llvm" \
&& clang --version \
&& ln -sf "/opt/llvm/bin/clang" "/usr/local/bin/cc" \
&& ln -sf "/opt/llvm/bin/clang++" "/usr/local/bin/c++" \
&& /install-openblas "x86_64-unknown-linux-gnu" "/opt/openblas" \
&& /install-compression-libs "x86_64-unknown-linux-gnu" "/usr/local" \
&& rm /install-openblas \
&& printf 'int main(void) { return 0; }\n' | cc -x c -o /tmp/cc-check - \
&& rm /tmp/cc-check

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
&& for dir in /tmp/lints/*; do (cd "${dir}" && rustup toolchain install); done \
&& rustup toolchain install "$(cat /tmp/hawk-toolchain)" --profile minimal --component rustc-dev,llvm-tools-preview,rust-src \
&& rustup toolchain install "$(cat /tmp/wasm-threads-toolchain)" --profile minimal --component rust-src,clippy \
&& chmod -R a+w "${RUSTUP_HOME}" "${CARGO_HOME}" \
&& rm -rf /install-rust /tmp/rust /tmp/lints /tmp/hawk-toolchain /tmp/wasm-threads-toolchain

# mise installs every tool of .config/mise.toml at the URL and sha256 of .config/mise.lock and
# links their executables into /usr/local/bin. The cargo registry that source
# builds fetch is a BuildKit cache, so it stays out of the image.
COPY .config/mise.toml .config/mise.lock /tmp/mise/project/
RUN --mount=type=cache,target=/usr/local/cargo/registry,sharing=locked \
  set -euxo pipefail >/dev/null \
&& export MISE_TRUSTED_CONFIG_PATHS="/tmp/mise/project" \
&& mise -C "/tmp/mise/project" install \
&& for dir in $(mise -C "/tmp/mise/project" bin-paths); do \
  find -L "${dir}" -mindepth 1 -maxdepth 1 -type f -executable -exec ln -sf -t "/usr/local/bin" {} + ; \
done \
&& chmod -R a+rX "${MISE_DATA_DIR}" \
&& rm -rf "/tmp/mise" \
&& just --version \
&& cargo nextest --version \
&& cargo dylint --version

# Emscripten builds the WebAssembly back end of the web2 app for
# wasm32-unknown-emscripten: rustc links through emcc, which finds Node on PATH.
# OpenBLAS comes prebuilt from binarylandia/build_openblas, and liblzma and
# libzstd from binarylandia/staticbourg, for this Emscripten version, because
# Emscripten gives no ABI guarantee between versions; an Emscripten upgrade
# needs new releases of both first. The single-threaded and the threaded
# OpenBLAS each get their own prefix, because wasm-ld rejects the
# single-threaded archive in a threaded link and the compression libraries put
# their prefix on every link line. Only the emcc front end goes on PATH, so the
# clang of Emscripten stays out of native builds. emcc runs the configuration
# .emscripten as Python at start: it names the LLVM and Binaryen of the archive
# and drops the include paths of the native GCC, which the clang of emcc would
# search too, mixing the glibc headers of the host into WebAssembly compiles.
# The cache holds the prebuilt system libraries, the threaded ones included,
# and stays writable for its lock and sanity files.
ENV PATH="/opt/emscripten/emscripten:${PATH}"
COPY dev/docker/files/install-openblas-emscripten /
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
&& /install-openblas-emscripten "st" "/usr/local/wasm32-unknown-emscripten-openblas" \
&& /install-openblas-emscripten "mt" "/usr/local/wasm32-unknown-emscripten-openblas-mt" \
&& /install-compression-libs "wasm32-unknown-emscripten" "/usr/local/wasm32-unknown-emscripten" \
&& chmod -R a+rwX "/opt/emscripten/emscripten/cache" \
&& rm -rf /install-openblas-emscripten /install-compression-libs /install-image-tools /mise-image "/tmp/mise"

# dev/docker/run runs as the host user with HOME=/tmp/home.
RUN set -euxo pipefail >/dev/null \
&& mkdir -p "/tmp/home" \
&& chmod 1777 "/tmp/home"
