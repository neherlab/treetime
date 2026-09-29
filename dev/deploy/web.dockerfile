# Runtime image of the web app: the static treetime-server binary, the built UI,
# and the example datasets. It only copies artifacts that dev/deploy/build-web-image
# builds beforehand, so building it takes seconds and a local build equals the CI build.
FROM gcr.io/distroless/static-debian12:nonroot@sha256:afa5c872c891853ca7fcf1f12c3edb23f7eeef36189728842dd51042ff57f7ab

ARG TARGET=x86_64-unknown-linux-musl

COPY .out/treetime-server-${TARGET} /app/treetime-server
COPY packages/app-web/dist /app/dist
COPY data /app/data

# /app mirrors the repository root: the example configs name their inputs
# relative to it (data/<dataset>/...), and the server resolves a relative input
# that is not inside the data directory from its working directory.
WORKDIR /app

ENV HOST=0.0.0.0
ENV PORT=3100
ENV STATIC_DIR=/app/dist

# The nonroot user of the distroless image.
USER 65532:65532
EXPOSE 3100
ENTRYPOINT ["/app/treetime-server"]
