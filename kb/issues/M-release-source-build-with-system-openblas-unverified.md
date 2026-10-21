# Source build with a system OpenBLAS unverified

## Summary

The shipped builds link a static OpenBLAS and the flang runtime from binarylandia/build_openblas inside the project containers. Whether a plain `cargo build --release` outside the containers links against a system OpenBLAS (found by `pkg-config` or the default library path) is not verified, documented, or tested. Packagers (bioconda, conda-forge, Debian, Nix, Homebrew) build from source with their own OpenBLAS, so the package channels depend on this path.

## Evidence

- `ndarray-linalg` uses the `openblas-system` feature ([Cargo.toml#L118](../../Cargo.toml#L118)) through the `numrs/openblas-src` fork ([Cargo.toml#L107](../../Cargo.toml#L107)).
- The containers replace the `openblas` link with the static archives (`dev/docker/files/install-openblas`) and add the archive to the target rustflags (`.cargo/config.toml`); no build runs without this override.

> [!IMPORTANT]
> **Investigation required.** Build the CLI on a clean Debian or Ubuntu with `libopenblas-dev` and on macOS with Homebrew OpenBLAS, without the containers, and record what fails (link names, LAPACK symbols, the Fortran runtime).

## Impact

- "Anyone can build it" holds only inside the project containers.
- Package recipes cannot be written until the plain build works.

## Potential solutions

- Make the plain build work with a system OpenBLAS, document it in the developer guide, and run it in CI on one Linux distribution and on macOS.
