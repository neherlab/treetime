# Cargo does not enforce the minimum publish age

## Summary

`.cargo/config.toml` sets `registry.global-min-publish-age = "7 days"` so that the resolver refuses crate versions published less than seven days ago. The pinned toolchain in `rust-toolchain.toml` is Rust 1.95.0, and its Cargo does not know this setting. `cargo update`, `cargo add`, and `cargo generate-lockfile` therefore resolve to the newest compatible versions, including versions published minutes earlier. The seven-day policy holds for Rust dependencies only when a developer checks publish dates by hand.

The JavaScript side is enforced: Bun reads `minimumReleaseAge` in `bunfig.toml` and refuses newer packages.

## Evidence

- Cargo added the setting as the unstable `-Zmin-publish-age` flag in [rust-lang/cargo#17012](https://github.com/rust-lang/cargo/pull/17012), merged 2026-06-18, and stabilized it for Rust 1.100 in [rust-lang/cargo#17335](https://github.com/rust-lang/cargo/pull/17335). Rust 1.95 contains neither.
- The Dylint toolchain (`nightly-2026-05-28`) is also older than the unstable flag, so Dylint runs cannot enforce the age either. Nightly Cargo reads the `[unstable]` table and does not know the key, so every Dylint build prints `warning: unused config key 'unstable.min-publish-age'`. The Dylint nightly is pinned to the toolchain of the `clippy_utils` revision that Dylint 6.0.4 is built against, so moving it past the flag means upgrading Dylint and every lint workspace together.
- Lockfile refreshes in this repository locked versions younger than seven days, for example `libredox 0.1.25` in the vendored Dylint workspace. Each one had to be held back by hand with `cargo update --precise`.

## Impact

- A compromised crate release can enter `Cargo.lock` on the next lockfile update, which defeats the purpose of the delay.
- Rust and JavaScript dependencies follow different rules: a too-new Bun package fails the install, while a too-new crate is accepted without a warning.

## Potential solutions

- Raise `rust-toolchain.toml` to Rust 1.100 or later once it is released. The resolver then honors `registry.global-min-publish-age` with no other change. This covers new resolutions only: Cargo does not check versions already in `Cargo.lock` ([rust-lang/cargo#17246](https://github.com/rust-lang/cargo/issues/17246)).
- Add a read-only check to `just check-all` that compares every registry entry in `Cargo.lock` with its crates.io publish time and fails on versions younger than seven days. This also covers lockfiles written by older Cargo versions or edited by hand, but it needs network access, so it fits `just audit` better than the offline gates.
