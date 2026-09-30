# kache does not seed new checkouts under the `.build/` target layout

## Summary

kache 0.27 and later copy the registry and git dependency units of a new checkout from another checkout's target directory, so Cargo compiles only the workspace crates ("target seeding", on by default). Seeding never runs for TreeTime, so every new worktree restores each dependency from the cache one compiler call at a time. Builds stay correct, but a new worktree does more work than it needs.

## Cause

kache derives a target's workspace from the target directory's parent and accepts it only when that parent holds a `Cargo.toml` (`verified_workspace_root()` in kache `src/args.rs`). TreeTime's target directories are `.build/host/cargo` and `.build/container/cargo` (`CARGO_TARGET_DIR` in the `justfile`), whose parents are `.build/host` and `.build/container`. kache then falls back to the working directory of a compiled package. `kache targets` shows both test targets of a two-checkout trial recorded under the registry source of `tikv-jemalloc-sys` instead of their checkouts. Seeding looks up the donor by the new checkout's `Cargo.lock`, and it finds none under that path.

The same misattribution affects kache's removal of orphaned targets. A target counts as orphaned once its recorded workspace path is gone (`workspace_is_gone()` in kache `src/cli.rs`). Cargo's cache cleanup can remove an unused registry source while the checkout still exists. `.cargo/config.toml` turns that removal off (`KACHE_AUTO_CLEAN_ORPHANED_TARGETS`), because every target directory lives under its checkout's `.build/` and goes away with it.

## Evidence

- Trial with kache 0.28.1 in the dev container, two fresh copies of the tree, `just test-rs --no-run` with a new store:
  - First copy: 173 s. Second copy: 51 s
  - The second copy's `kache report --last-build` lists 549 of 551 crates as cache hits. Seeded units would not reach the wrapper
- kache v0.27.0 release notes describe seeding: host `debug` profile only, needs a running daemon ([kache v0.27.0](https://github.com/kunobi-ninja/kache/releases/tag/v0.27.0))

## Resolution options

- Place each target directory directly under the checkout root, for example `CARGO_TARGET_DIR=<checkout>/.build-container`, so its parent is the workspace. This changes the build layout that `dev/docker/run`, the justfile, and the lint target directories share
- Ask kache upstream to take the workspace from Cargo's environment (`CARGO_MANIFEST_DIR` of a workspace member, or the `Cargo.lock` location) when the target directory is external
