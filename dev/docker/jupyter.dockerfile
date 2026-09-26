# syntax=docker/dockerfile:1
# check=experimental=all
FROM debian:12.14-slim@sha256:60eac759739651111db372c07be67863818726f754804b8707c90979bda511df

SHELL ["bash", "-euxo", "pipefail", "-c"]

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
  lsb-release \
  parallel \
  pigz \
  pixz \
  pkg-config \
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

COPY --link "dev/docker/files/fetch" "dev/docker/files/checksums" "dev/docker/files/create-user" "/"
RUN /create-user

# The Python environment comes from pixi.toml and pixi.lock, the same files that
# host users run with `pixi`. The pixi version follows .config/mise.toml.
ENV PIXI_PROJECT="/opt/python"
ARG PIXI_VERSION="0.81.0"
COPY --link "pixi.toml" "pixi.lock" "${PIXI_PROJECT}/"
RUN set -euxo pipefail >/dev/null \
&& /fetch "https://github.com/prefix-dev/pixi/releases/download/v${PIXI_VERSION}/pixi-x86_64-unknown-linux-musl" "/usr/local/bin/pixi" \
&& chmod +x "/usr/local/bin/pixi" \
&& pixi install --locked --manifest-path "${PIXI_PROJECT}/pixi.toml" \
&& pixi clean cache --yes \
&& chmod -R a+rX "${PIXI_PROJECT}"

COPY --link "dev/docker/files/.jupyter/" "/.jupyter"
COPY --link "dev/docker/files/start-jupyter" "/"
COPY --link --chmod=0755 "dev/docker/files/usr/bin/treetime" "/usr/bin/treetime"

USER ${USER}

ENV PATH="${PIXI_PROJECT}/.pixi/envs/default/bin:${PATH}"
ENV OPENBLAS_NUM_THREADS="1"
ENV XDG_CACHE_HOME="${HOME}/.cache/"
ENV MPLBACKEND="Agg"

# Import matplotlib the first time to build the font cache.
RUN set -euxo pipefail >/dev/null \
&& python -c "import matplotlib.pyplot"
