# Cargo does not enforce the minimum publish age on its own

## Summary

Cargo can refuse crate versions published less than seven days ago through `registry.global-min-publish-age`. The pinned toolchain, Rust 1.98.1, has this only as the unstable `-Zmin-publish-age` feature; it becomes stable in Rust 1.100, and a stable Cargo warns on every command when the key is set in `.cargo/config.toml`. A plain `cargo update`, `cargo add`, or `cargo generate-lockfile` therefore resolves to the newest compatible versions, including versions published minutes earlier.

The JavaScript side is enforced: Bun reads `minimumReleaseAge` in `bunfig.toml` and refuses newer packages.

## Current mitigation

- `just deps-update` and `just deps-upgrade` run Cargo with `RUSTC_BOOTSTRAP=1 -Zmin-publish-age` and the setting passed with `--config`, which enables the unstable feature on the stable toolchain for the resolution only
- `just deps-age` (`dev/crate-age`, part of `just audit`) reads every crates.io entry of `Cargo.lock` and of the lint workspace lockfiles, and fails on a version younger than seven days. It needs network access, so it is not part of the offline gates
- Builds run with `--locked`, so no other recipe changes `Cargo.lock`

## Remaining gap

- A direct `cargo update` or `cargo add` outside the recipes still accepts a young release. `just deps-age` finds it, but only when someone runs it
- `cargo install` of the `cargo:` tools in `mise.toml` does not apply the age check at all, in any Cargo version
- Git dependencies in `[patch.crates-io]` have no publish time

## Resolution

Raise `rust-toolchain.toml` to Rust 1.100 or later once it is at least seven days old, set `registry.global-min-publish-age = "7 days"` in `.cargo/config.toml`, and drop the flags of the dependency recipes (`cargo_min_age` in the justfile). Cargo does not check versions already in `Cargo.lock` ([rust-lang/cargo#17246](https://github.com/rust-lang/cargo/issues/17246)), so `just deps-age` keeps its value for lockfiles written by other tools.

## Evidence

- Cargo added the setting as the unstable `-Zmin-publish-age` flag in [rust-lang/cargo#17012](https://github.com/rust-lang/cargo/pull/17012) and stabilized it for Rust 1.100 in [rust-lang/cargo#17335](https://github.com/rust-lang/cargo/pull/17335)
- Cargo treats `RUSTC_BOOTSTRAP=1` as a development channel that allows unstable features (`src/cargo/core/features.rs`, `channel()`)
