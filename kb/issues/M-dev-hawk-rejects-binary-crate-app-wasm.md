# `just hawk` stops on the binary crate `app-wasm`

## Summary

`just hawk` stops before any check with `unknown --exclude-crate value(s): app_wasm`, so the check of unnecessary public items (cargo-hawk) does not run.

## Cause

The recipe passes every crate of the `public_api_crates` variable of the `justfile` as `--exclude-crate=<crate>`. The list includes `app_wasm`, but `app-wasm` is a binary crate: it has only `packages/app-wasm/src/main.rs` and no library target. `cargo hawk --exclude-crate` accepts only library crates and rejects the whole command line for an unknown value.

The recipe `dylint` passes the same list to `pub-unused-report` as `--exclude-crate`.

## Resolution

- Remove `app_wasm` from `public_api_crates`: a binary crate has no public API that other crates use, so cargo-hawk has nothing to exclude there. Check that `pub-unused-report` of `just dylint` still accepts the list, then run `just hawk`
