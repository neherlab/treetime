# syntax=docker/dockerfile:1
# check=experimental=all
#
# Production image neherlab/treetime without a base: only the static musl build
# of the CLI as /usr/bin/treetime, which reads and writes local files only.
# dev/docker/prod-images stages the binary of each platform under
# .out/docker/<arch>/musl/.
FROM scratch

ARG TARGETARCH
COPY ".out/docker/${TARGETARCH}/musl/treetime" "/usr/bin/treetime"

ENV PATH="/usr/bin"
