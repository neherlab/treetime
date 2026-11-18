# `just dylint` cannot build its driver in the container: pkg-config finds no OpenSSL

## Summary

In a checkout without a built dylint driver (`.cache/cargo/dylint_drivers/` empty), `just dylint` stops before any lint runs. `cargo dylint` builds the driver for the toolchain of the lint libraries (`nightly-2026-05-28`) in a temporary directory, and the build of `openssl-sys` (a dependency through `curl` and `libgit2-sys`) fails with `Could not find openssl via pkg-config`. The recipe passes `--quiet`, so its log shows only `command failed: ... "cargo" "build"`, followed by one `pub-unused-report: no record for ...` line per target and `the report is incomplete and nothing is judged`.

## Cause

The image writes `PKG_CONFIG_PATH_x86_64_unknown_linux_gnu` and `PKG_CONFIG_LIBDIR_x86_64_unknown_linux_gnu` into the `[env]` table of `/.cargo/config.toml` (`fn pkg_config_env()` in `dev/docker/files/install-image-libraries`). Both name only the pkg-config folders of the prebuilt OpenBLAS, liblzma and libzstd archives, so the static builds of the project never link a system library. Cargo applies that table to every build below `/`, so the driver build in `/tmp` gets the same paths, and pkg-config does not search `/usr/lib/x86_64-linux-gnu/pkgconfig`, where the image has `openssl.pc` (`libssl-dev` is installed).

## Fix direction

Let the driver build find the system OpenSSL without widening the library search of the project builds, for example by setting `OPENSSL_DIR` or `OPENSSL_LIB_DIR` and `OPENSSL_INCLUDE_DIR` in the `dylint` recipe or the image, or by building the driver with an environment without the `[env]` paths of the image.

## Validation

- `just dylint` in a new checkout builds the driver and reports lint findings
