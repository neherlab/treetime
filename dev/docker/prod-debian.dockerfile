# syntax=docker/dockerfile:1
# check=experimental=all
#
# Production image neherlab/treetime on Debian: the glibc build of the CLI as
# /usr/bin/treetime, with bash, ca-certificates, curl, procps, and wget, because
# workflow managers such as Nextflow need bash and ps in task containers.
# dev/docker/prod-images builds it on the Debian base of
# dev/docker/prod-bases.json, passed as BASE_IMAGE with its digest, and stages
# the binary of each platform under .out/docker/<arch>/gnu/.
ARG BASE_IMAGE
FROM ${BASE_IMAGE}
SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  bash \
  ca-certificates \
  curl \
  procps \
  wget \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

ARG TARGETARCH
COPY ".out/docker/${TARGETARCH}/gnu/treetime" "/usr/bin/treetime"
