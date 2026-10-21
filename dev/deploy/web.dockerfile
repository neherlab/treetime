# Runtime image of the web app: the static treetime-server binary, the built UI,
# and the example datasets. It only copies artifacts that dev/deploy/build-web-image
# builds beforehand, so building it takes seconds and a local build equals the CI build.

# CA bundle from the current Debian package; the one of the distroless image is older.
FROM debian:12.15-slim@sha256:7ff6ba31527b9e336162fa061746c366c4b2032f4d2cff2078bea43f78efe3f8 AS certs

SHELL ["bash", "-euxo", "pipefail", "-c"]

RUN set -euxo pipefail >/dev/null \
&& export DEBIAN_FRONTEND=noninteractive \
&& apt-get update -qq \
&& apt-get install --no-install-recommends --yes -qq \
  ca-certificates \
>/dev/null \
&& apt-get clean autoclean >/dev/null \
&& apt-get autoremove --yes >/dev/null \
&& rm -rf /var/lib/apt/lists/* /var/cache/apt/archives/*

FROM gcr.io/distroless/static-debian12:nonroot@sha256:afa5c872c891853ca7fcf1f12c3edb23f7eeef36189728842dd51042ff57f7ab

ARG TARGET=x86_64-unknown-linux-musl

COPY --from=certs /etc/ssl/certs/ca-certificates.crt /etc/ssl/certs/ca-certificates.crt
COPY .out/treetime-server-${TARGET} /app/treetime-server
COPY packages/app-web/dist /app/dist
COPY data /app/data

# The server resolves a relative input path inside the examples folder first,
# then from this working directory.
WORKDIR /app

ENV HOST=0.0.0.0
ENV PORT=3100
ENV STATIC_DIR=/app/dist

# The nonroot user of the distroless image.
USER 65532:65532
EXPOSE 3100
ENTRYPOINT ["/app/treetime-server"]
