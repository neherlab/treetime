# syntax=docker/dockerfile:1
# check=experimental=all
#
# Production image neherlab/treetime on Alpine: the static musl build of the
# CLI as /usr/bin/treetime. dev/docker/prod-images builds it on the Alpine base
# of dev/docker/prod-bases.json, passed as BASE_IMAGE with its digest, and
# stages the binary of each platform under .out/docker/<arch>/musl/.
ARG BASE_IMAGE
FROM ${BASE_IMAGE}

ARG TARGETARCH
COPY ".out/docker/${TARGETARCH}/musl/treetime" "/usr/bin/treetime"
