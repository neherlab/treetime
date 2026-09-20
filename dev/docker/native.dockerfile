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
  libfontconfig-dev \
  libssl-dev \
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

COPY --link "dev/docker/files/install-seqkit" "/"
RUN /install-seqkit

COPY --link "dev/docker/files/install-iqtree" "/"
RUN /install-iqtree

COPY --link "dev/docker/files/install-kache" "/"
RUN /install-kache

COPY --link "dev/docker/files/install-nodejs" "/"
RUN /install-nodejs

COPY --link "dev/docker/files/install-electron-deps" "/"
RUN /install-electron-deps

ENV GTK_THEME="Adwaita:dark"

RUN set -euxo pipefail >/dev/null \
&& mkdir -p "/etc/gtk-3.0" \
&& printf '[Settings]\ngtk-application-prefer-dark-theme=1\n' > "/etc/gtk-3.0/settings.ini"

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

COPY --link "dev/docker/files/install-dylint" "/"
COPY --link --chown="${UID}:${GID}" "dev/lints/dylint/rust-toolchain.toml" "/tmp/lints/dylint/rust-toolchain.toml"
RUN set -euxo pipefail >/dev/null \
&& /install-dylint \
&& rm -rf "/tmp/lints"

COPY --link "dev/docker/files/install-hawk" "dev/docker/files/hawk-toolchain" "/"
RUN set -euxo pipefail >/dev/null \
&& /install-hawk


# Developer tooling via mise, pinned and checksummed in mise.lock. This layer is
# last so the expensive toolchain layers above stay cached. mise.toml / mise.lock
# are inputs to dev/docker/checksum, so a tool change rebuilds the image. At run
# time the mise shims resolve the same versions from the global config below.
#
# mise fetches `ubi:` and `github:` tools through the GitHub API, which caps
# anonymous callers at 60 requests/hour/IP and fails the whole build once that is
# spent. dev/docker/run passes a token through a BuildKit secret so mise
# authenticates; `env=GITHUB_TOKEN` exposes it only for this RUN, so it never
# reaches an image layer or the build history. The secret is optional: with no
# token GITHUB_TOKEN stays unset and mise falls back to anonymous access.
ENV MISE_DATA_DIR="${HOME}/.local/share/mise"
ENV MISE_GLOBAL_CONFIG_FILE="${HOME}/tools/mise.toml"
ENV MISE_TRUSTED_CONFIG_PATHS="/workdir:${HOME}/tools"
ENV MISE_NOT_FOUND_AUTO_INSTALL="0"
ENV PATH="${HOME}/.local/bin:${MISE_DATA_DIR}/shims:${PATH}"
ARG MISE_VERSION="v2026.9.11"
COPY --link --chown="${UID}:${GID}" "mise.toml" "mise.lock" "${HOME}/tools/"
RUN --mount=type=secret,id=github_token,env=GITHUB_TOKEN set -euxo pipefail >/dev/null \
&& curl -fsSL https://mise.run | MISE_INSTALL_PATH="${HOME}/.local/bin/mise" MISE_VERSION="${MISE_VERSION}" sh \
&& mise trust "${HOME}/tools/mise.toml" \
&& mise install \
&& mise reshim \
&& mise ls \
&& just --version
