# CDS selection is read from the settings in two places

Two functions read the list of CDS names from the `cds_selection` key of the run settings, with the same code and each with its own copy of the key:

- [packages/app-commands/src/commands/ancestral/aa_node_data.rs#L16](../../packages/app-commands/src/commands/ancestral/aa_node_data.rs#L16) and [packages/app-commands/src/commands/ancestral/aa_node_data.rs#L32-L39](../../packages/app-commands/src/commands/ancestral/aa_node_data.rs#L32-L39): `CDS_SELECTION_KEY` and `pub fn cds_selection_names()`
- [packages/app-commands/src/commands/shared/run_paths.rs#L19](../../packages/app-commands/src/commands/shared/run_paths.rs#L19) and [packages/app-commands/src/commands/shared/run_paths.rs#L251-L258](../../packages/app-commands/src/commands/shared/run_paths.rs#L251-L258): `CDS_SELECTION_KEY` and `fn listed_cds_names()`

When the key or the shape of the setting changes, a change of one copy leaves the run paths and the amino-acid node data reading different CDS lists.

## Fix

Keep one reader of the CDS selection, in the module that owns the run settings, and call it from both places.

## Validation

`just test-rs -E 'package(app-commands)'` passes, and `rg 'CDS_SELECTION_KEY' packages` finds one definition.

## Proposed implementation plan

Both copies are byte-identical and the cited lines are current. A third caller, [packages/app-cli/src/cli/pipeline/safety.rs#L144](../../packages/app-cli/src/cli/pipeline/safety.rs#L144), already calls `cds_selection_names()`.

1. Keep `pub fn cds_selection_names()` and `CDS_SELECTION_KEY` in `packages/app-commands/src/commands/ancestral/aa_node_data.rs`, because the ancestral command owns the `cds_selection` setting and `run_paths.rs` already imports `cds_output_paths` and `translation_path` from that module
2. In `packages/app-commands/src/commands/shared/run_paths.rs`: delete `const CDS_SELECTION_KEY` and `fn listed_cds_names()`, add `cds_selection_names` to the existing `use crate::commands::ancestral::aa_node_data::{...}` import, and call it in `fn check_command_paths()`
3. Tests: no new test. The behavior is unchanged, and the existing cases of `test_aa_node_data.rs` (`listed_names_first`, `empty_selection`) and the `cds_selection` cases of `packages/app-cli/src/cli/pipeline/safety.rs` and `packages/app-server/src/__tests__/test_confine.rs` cover both callers

### Plan validation

- `rg 'CDS_SELECTION_KEY' packages` finds one definition and `rg 'listed_cds_names' packages` finds nothing
- `just lint-rs` and `just test-rs` pass
