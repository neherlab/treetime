# syntax=docker/dockerfile:1
# check=experimental=all
#
# Development image: Rust and Bun toolchains, lint and test tools, and the
# libraries the desktop app needs. dev/docker/run builds it and runs every
# command as the host user, with HOME=/tmp/home and the cargo home in the
# checkout. Ubuntu 24.04 keeps the glibc of the binaries built here no newer than
# on current hosts, so profiling builds run on the host.
FROM ubuntu:noble-20260905@sha256:a053cbffda9d424679c103c5b4f452297efc3774a1e491c110289f726fbb5d34
SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  ca-certificates \
  curl \
  file \
  git \
  libc6-dev \
  libfontconfig-dev \
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

# GCC 14 is the C compiler and link driver; it also provides the static
# libgfortran that the OpenBLAS archive needs. It does not know the Debian
# multiarch layout, so the include and library paths name it.
# LLVM provides libclang for bindgen.
ENV C_INCLUDE_PATH="/usr/include/x86_64-linux-gnu"
ENV CPLUS_INCLUDE_PATH="/usr/include/x86_64-linux-gnu"
ENV LIBRARY_PATH="/usr/local/lib64:/usr/local/lib:/usr/lib/x86_64-linux-gnu"
ENV LD_LIBRARY_PATH="/usr/local/lib64:/usr/local/lib"
COPY dev/docker/files/install-gcc dev/docker/files/install-llvm dev/docker/files/install-openblas /
RUN set -euxo pipefail >/dev/null \
&& /install-gcc "/usr/local" \
&& /install-llvm "/usr/local" \
&& /install-openblas "x86_64-unknown-linux-gnu" "/usr/local" \
&& rm /install-gcc /install-llvm /install-openblas

# The Rust toolchain of rust-toolchain.toml, then the pinned nightly of the lint
# libraries and the toolchain that cargo-hawk is built against.
ENV RUSTUP_HOME="/usr/local/rustup"
ENV CARGO_HOME="/usr/local/cargo"
ENV PATH="/usr/local/cargo/bin:${PATH}"
COPY dev/docker/files/install-rust /
COPY rust-toolchain.toml /tmp/rust/
COPY dev/lints/dylint-custom/rust-toolchain.toml /tmp/lints/dylint-custom/
COPY dev/lints/dylint-mordant/rust-toolchain.toml /tmp/lints/dylint-mordant/
COPY dev/lints/dylint-trailofbits/rust-toolchain.toml /tmp/lints/dylint-trailofbits/
COPY dev/docker/files/hawk-toolchain /tmp/
# hadolint ignore=DL3003
RUN set -euxo pipefail >/dev/null \
&& /install-rust "/tmp/rust" \
&& for dir in /tmp/lints/*; do (cd "${dir}" && rustup toolchain install); done \
&& rustup toolchain install "$(cat /tmp/hawk-toolchain)" --profile minimal --component rustc-dev,llvm-tools-preview,rust-src \
&& chmod -R a+w "${RUSTUP_HOME}" "${CARGO_HOME}" \
&& rm -rf /install-rust /tmp/rust /tmp/lints /tmp/hawk-toolchain

# mise installs every tool of mise.toml at the URL and sha256 of mise.lock and
# links their executables into /usr/local/bin. The cargo registry that source
# builds fetch is a BuildKit cache, so it stays out of the image.
ENV MISE_DATA_DIR="/opt/mise"
ENV MISE_CACHE_DIR="/tmp/mise/cache"
ENV MISE_STATE_DIR="/tmp/mise/state"
COPY mise.toml mise.lock /tmp/mise/project/
RUN --mount=type=cache,target=/usr/local/cargo/registry,sharing=locked \
  set -euxo pipefail >/dev/null \
&& /fetch "https://github.com/jdx/mise/releases/download/v2026.9.10/mise-v2026.9.10-linux-x64-musl.tar.gz" "/tmp/mise.tar.gz" \
&& tar -xzf "/tmp/mise.tar.gz" --strip-components=2 -C "/usr/local/bin" "mise/bin/mise" \
&& export MISE_TRUSTED_CONFIG_PATHS="/tmp/mise/project" \
&& mise -C "/tmp/mise/project" install \
&& for dir in $(mise -C "/tmp/mise/project" bin-paths); do \
  find -L "${dir}" -mindepth 1 -maxdepth 1 -type f -executable -exec ln -sf -t "/usr/local/bin" {} + ; \
done \
&& chmod -R a+rX "${MISE_DATA_DIR}" \
&& rm -rf "/tmp/mise" "/tmp/mise.tar.gz" \
&& just --version \
&& cargo nextest --version \
&& cargo dylint --version

# dev/docker/run runs as the host user with HOME=/tmp/home.
RUN set -euxo pipefail >/dev/null \
&& mkdir -p "/tmp/home" \
&& chmod 1777 "/tmp/home"
