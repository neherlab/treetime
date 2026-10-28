# syntax=docker/dockerfile:1
# check=experimental=all
#
# Cross image of the WebAssembly back end of the web2 app
# (wasm32-unknown-emscripten), for `just build-web2` in the web deploy:
# Emscripten, the Rust toolchain of rust-toolchain.toml with the nightly of the
# threaded module, the prebuilt WebAssembly libraries, and the tools that the
# build runs. It leaves out the native toolchains, the lint toolchains, and the
# development tools of the native image, so it builds, pushes, and pulls in a
# fraction of the time. The GCC of Debian compiles and links the build scripts
# and procedural macros, which only run during the build.
#
# The layers go from the inputs that change least to those that change most:
# the Rust toolchains (rust-toolchain.toml), Emscripten and the libraries
# (.config/mise.image.lock), then the tools (.config/mise.lock).
FROM debian:13.7-slim@sha256:62bd423aeee912639ca8ee681c2b762cb3e25137e0a92b20be455898a9f034ad

SHELL ["bash", "-euxo", "pipefail", "-c"]

# python3 runs emcc.
RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  bzip2 \
  ca-certificates \
  curl \
  file \
  gcc \
  git \
  libc6-dev \
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

COPY dev/docker/files/fetch dev/docker/files/checksums /

# The Rust toolchain of rust-toolchain.toml, which includes the
# wasm32-unknown-emscripten target, and the nightly of the threaded module,
# which rebuilds std with atomics from rust-src (-Zbuild-std).
ENV RUSTUP_HOME="/usr/local/rustup"
ENV CARGO_HOME="/usr/local/cargo"
ENV PATH="/usr/local/cargo/bin:${PATH}"
COPY dev/docker/files/install-rust /
COPY rust-toolchain.toml /tmp/rust/
COPY dev/docker/files/wasm-threads-toolchain /tmp/
RUN set -euxo pipefail >/dev/null \
&& /install-rust "/tmp/rust" \
&& rustup toolchain install "$(cat /tmp/wasm-threads-toolchain)" --profile minimal --component rust-src,clippy \
&& chmod -R a+w "${RUSTUP_HOME}" "${CARGO_HOME}" \
&& rm -rf /install-rust /tmp/rust /tmp/wasm-threads-toolchain

# mise installs Emscripten and the libraries of .config/mise.image.toml here and
# the tools of .config/mise.toml further down, each at the URLs and checksums of
# its lockfile.
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

# Emscripten, as in dev/docker/native.dockerfile: rustc links through emcc,
# which finds Node on PATH. The configuration .emscripten names the LLVM and
# Binaryen of the archive. The cache holds the prebuilt system libraries, the
# threaded ones included, and stays writable for its lock and sanity files.
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
&& chmod -R a+rwX "/opt/emscripten/emscripten/cache"

# OpenBLAS, liblzma, and libzstd of the single-threaded and the threaded module:
# the prebuilt archives of .config/mise.image.toml, which cargo links through
# pkg-config with the [env] table that the install writes to
# /.cargo/config.toml and /.cargo/wasm-mt-env.toml. The install fails when an
# archive was built with another Emscripten version than this image has.
COPY dev/docker/files/install-image-libraries dev/lib/dist-flags.sh /
RUN set -euxo pipefail >/dev/null \
&& /install-image-libraries "wasm32-unknown-emscripten" \
&& rm -rf /install-image-libraries /dist-flags.sh /install-image-tools /mise-image "/tmp/mise"

# The tools of `just build-web2`, at the versions, URLs, and sha256 of
# .config/mise.lock: just, which runs every recipe through `mise exec`; bun and
# node for the Vite build and emcc; wasm-bindgen, which emcc runs; protoc, which
# prost-build runs for the UShER protobuf schemas; and kache, the optional
# compiler cache. The system configuration of mise disables the other tools of
# .config/mise.toml, so `mise exec` neither warns about them nor installs them.
COPY .config/mise.toml .config/mise.lock /tmp/mise/project/
RUN set -euxo pipefail >/dev/null \
&& export MISE_TRUSTED_CONFIG_PATHS="/tmp/mise/project" \
&& set -- just bun node wasm-bindgen protoc "github:kunobi-ninja/kache" \
&& mise -C "/tmp/mise/project" install "$@" \
&& mise -C "/tmp/mise/project" bin-paths "$@" | while IFS= read -r dir; do \
  find -L "${dir}" -mindepth 1 -maxdepth 1 -type f -executable -exec ln -sf -t "/usr/local/bin" {} + || exit 1; \
done \
&& disabled="$(mise -C "/tmp/mise/project" ls --missing --json | python3 -c 'import json, sys; print(json.dumps(sorted(json.load(sys.stdin))))')" \
&& mkdir -p "/etc/mise" \
&& printf '[settings]\nexec_auto_install = false\ndisable_tools = %s\n' "${disabled}" >"/etc/mise/config.toml" \
&& chmod -R a+rX "${MISE_DATA_DIR}" \
&& rm -rf "/tmp/mise" \
&& cat "/etc/mise/config.toml" \
&& just --version \
&& bun --version \
&& node --version \
&& wasm-bindgen --version \
&& protoc --version \
&& kache --version \
&& emcc --version

# dev/docker/run runs as the host user with HOME=/tmp/home.
RUN set -euxo pipefail >/dev/null \
&& mkdir -p "/tmp/home" \
&& chmod 1777 "/tmp/home"
