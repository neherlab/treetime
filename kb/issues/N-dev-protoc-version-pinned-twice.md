# protoc version is pinned in two places

## Summary

The protobuf compiler, which `prost-build` runs for `packages/util-usher-mat`, has two version owners. `.config/mise.toml` pins `protoc` for the host and the development container. `dev/docker/files/install-protobuf` hardcodes the same version for the cross-compilation images (`dev/docker/cross-linux.dockerfile`, `dev/docker/cross-darwin.dockerfile`). An upgrade that changes one file leaves the other on the old version, so native builds and cross builds generate code with different compilers.

## Evidence

- `.config/mise.toml`: `protoc = "28.3"`
- `dev/docker/files/install-protobuf`: `version="28.3"`, verified by its own sha256 entry in `dev/docker/files/checksums`, separate from the one in `.config/mise.lock`

## Resolution options

- Install protoc in the cross images through mise from `.config/mise.toml` and `.config/mise.lock`, as the development container does, and delete `install-protobuf`
- Read the version in `install-protobuf` from `.config/mise.toml` and verify the archive against the checksum in `.config/mise.lock`
