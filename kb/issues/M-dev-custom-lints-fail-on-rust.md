# Custom lints fail on the `rust` branch

## Summary

`just dylint` fails on `rust`, and so does the CI job "Check (custom lints)" from commit `a013d8997` on. The lint libraries of `dev/lints/` (`custom`, `mordant`, `trailofbits`) report 55 findings in the crates they reach, and five application crates are not linted at all, so their findings are unknown.

## Crates without a result

The `no_comments` lint of the `custom` library rejects the `///` lines in the `macro_rules!` body of `packages/app-commands/src/commands/shared/ambiguity.rs` (18 errors). The macro generates a `clap::Args` struct, but the lint sees the doc comments in the macro body, where no derive renders them. The errors stop the compilation of `app-commands` under the lint, so every crate that depends on it gets no lint run: `app-backend`, `app-cli`, `app-napi`, `app-server` and `app-wasm`. The workspace report of `pub_unused_in_workspace` then stops with `pub-unused-report: no record for lib app_backend ...` (one line per target) and `the report is incomplete and nothing is judged`.

`packages/app-commands/src/commands/shared/output_args.rs` has the same macro pattern and passes its help text as `#[doc = "..."]` attributes, which the lint accepts.

## Findings by category

- **Doc comments in a macro body** (`no_comments`, errors): `packages/app-commands/src/commands/shared/ambiguity.rs` (18)
- **Item order** (`topological_ordering`, items not ordered callers before callees): `packages/treetime-distribution/src/distribution_ops/convolve.rs`, `packages/treetime-distribution/src/distribution_ops/mass_domain.rs`, `packages/treetime-distribution/src/distribution_ops/mass_profile.rs`, `packages/treetime-graph/src/clades.rs`, `packages/treetime-simulation/src/coalescent.rs`, `packages/treetime-simulation/src/evolve.rs`, `packages/treetime/src/coalescent/knot_times.rs`, `packages/treetime/src/partition/marginal/discrete/partition.rs`, `packages/treetime/src/seq/sink.rs`, `packages/treetime/src/timetree/coalescent_timescale.rs`, `packages/treetime/src/timetree/inference/result.rs`, `packages/treetime/src/timetree/max_posterior_region.rs`
- **Construction that skips the checked constructor** (`mordant`, "built by hand here, skipping the check" and "assigned directly here"):
  - `HardApproachLaw`: `packages/treetime-grid/src/__tests__/grid_fn.rs` (5)
  - `NodePosterior`: `packages/treetime/src/coalescent/__tests__/test_mean_field.rs` (3), `packages/treetime/src/test_utils/timetree.rs` (2), `packages/treetime/src/timetree/__tests__/test_confidence_extract.rs`, `packages/treetime/src/timetree/__tests__/test_time_inference.rs`
  - `NodePosterior::time` assigned directly: `packages/treetime/src/timetree/__tests__/test_round.rs` (4), `packages/treetime/src/timetree/convergence/__tests__/test_node_times.rs`
- **Catch-all arm over an own enum** (`wildcard_over_own_enum`): `packages/treetime-distribution/src/distribution_core/distribution.rs` and `packages/treetime-distribution/src/distribution_ops/__tests__/test_mass_trim.rs` (`Distribution`), `packages/treetime-distribution/src/distribution_ops/mass_profile.rs` (`Cell`)
- **Repeated match** (`same_match_twice`, on `treetime_grid::Side`): `packages/treetime-distribution/src/distribution_ops/__tests__/test_divide.rs`
- **Discarded error** (`let Ok(..) = .. else` drops the error): `packages/treetime/src/partition/marginal/sparse/reconstruct.rs`
- **Renamed import**: `packages/treetime/src/timetree/max_posterior_region.rs`
- **Version suffix** (`versioned_name`): `test_augur_node_data_mugration_model_and_root_match_v0` in `packages/app-commands/src/commands/mugration/__tests__/test_augur_node_data.rs`. Here `v0` names the reference implementation of the test, TreeTime 0.x
- **Duplicate helper** (`mordant`, same signature and body): the test helper `function` in `packages/treetime-distribution/src/distribution_ops/__tests__/test_mass_profile_tails.rs` copies the one in `packages/treetime-distribution/src/distribution_ops/__tests__/test_mass_profile.rs`

## Resolution

- Pass the help text in `ambiguity.rs` as `#[doc = "..."]` attributes, as `output_args.rs` does, and check with `just gen` that the CLI reference and the schemas do not change
- Run `just dylint` again: the five crates are then linted, and their findings join the list above
- Fix each finding, or suppress one that detects no real defect with `#[cfg_attr(dylint_lib = "<library>", expect(<lint>, reason = ".."))]` on the narrowest item
