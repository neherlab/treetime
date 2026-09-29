# Two CI check groups fail on `rust`

The `cli.yml` run on `rust` at `400f4ab9` fails two check groups of `just check-all`. The nightly CLI release does not wait for these checks, but every pull request and push to `rust` reports red until they pass.

- **`Check (format and config lints)`**: `just shear` reports `async-stream` as unused, both in `packages/app-server/Cargo.toml` and in the workspace dependency table of the root `Cargo.toml`. No Rust source uses it. The failure reproduces locally with `just shear`
- **`Check (custom lints)`**: `just dylint` denies warnings from the custom lint libraries. Most findings are in the app crates (`app-commands`, `app-server`, `app-napi`), with some in `treetime-io` and `treetime`. The most frequent lints:
  - `topological_ordering`: items are not ordered callers-before-callees in a module
  - `try_io_result`: a function returns `std::io::Result`, which drops the file or path context
  - `result_defaulted`: `.ok()` replaces an error with a default and loses the cause
  - `spawn_handle_dropped`: the `JoinHandle` of a spawned task is dropped immediately
  - Also reported: `file_too_long`, `test_real_sleep`, `assert_in_loop`, `error_dropped_by_pattern`, `suggest_builder`, `unnamed_constant`, and a `from_value` call that copies the whole `Value` before parsing

`just dylint` is a slow recipe, so the dylint findings come from the CI log rather than from a local run.

## Resolution

Remove `async-stream` from both manifests. Fix the dylint findings, or suppress one with a justification where it does not point at a real defect.
