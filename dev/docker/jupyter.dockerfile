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

# pixi at the version, URL, and sha256 of .config/mise.lock, so the host and the
# image run the same pixi.
ENV MISE_DATA_DIR="/opt/mise"
ENV MISE_CACHE_DIR="/tmp/mise/cache"
ENV MISE_STATE_DIR="/tmp/mise/state"
COPY --link "dev/docker/files/install-mise" "dev/docker/files/mise-version" "/"
COPY --link ".config/mise.toml" ".config/mise.lock" "/tmp/mise/project/"
RUN set -euxo pipefail >/dev/null \
&& /install-mise "/usr/local/bin" \
&& export MISE_TRUSTED_CONFIG_PATHS="/tmp/mise/project" \
&& mise -C "/tmp/mise/project" install "pixi" \
&& ln -sf -t "/usr/local/bin" "$(mise -C "/tmp/mise/project" which pixi)" \
&& chmod -R a+rX "${MISE_DATA_DIR}" \
&& rm -rf "/tmp/mise" "/install-mise" "/mise-version" "/usr/local/bin/mise" \
&& pixi --version

# The Python environment comes from pixi.toml and pixi.lock, the same files that
# host users run with `pixi`.
ENV PIXI_PROJECT="/opt/python"
COPY --link "pixi.toml" "pixi.lock" "${PIXI_PROJECT}/"
RUN set -euxo pipefail >/dev/null \
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
