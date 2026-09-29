# syntax=docker/dockerfile:1
# check=experimental=all
FROM debian:12.14-slim@sha256:60eac759739651111db372c07be67863818726f754804b8707c90979bda511df

SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  bash \
  gifsicle \
  imagemagick \
  libxml2-utils \
  librsvg2-bin \
  nodejs \
  npm \
  potrace \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

ENV SVGO_DIR="/opt/svgo"
ENV PATH="${SVGO_DIR}/node_modules/.bin:${PATH}"
COPY --link "package.json" "package-lock.json" "${SVGO_DIR}/"
RUN set -euxo pipefail >/dev/null \
&& npm ci --prefix="${SVGO_DIR}" --ignore-scripts --omit=dev --no-audit --no-fund --cache="/tmp/npm" \
&& rm -rf "/tmp/npm"
