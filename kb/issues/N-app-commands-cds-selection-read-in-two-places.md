# CDS selection is read from the settings in two places

Two functions read the list of CDS names from the `cds_selection` key of the run settings, with the same code and each with its own copy of the key:

- [packages/app-commands/src/commands/ancestral/aa_node_data.rs#L16](../../packages/app-commands/src/commands/ancestral/aa_node_data.rs#L16) and [packages/app-commands/src/commands/ancestral/aa_node_data.rs#L27-L34](../../packages/app-commands/src/commands/ancestral/aa_node_data.rs#L27-L34): `CDS_SELECTION_KEY` and `pub fn cds_selection_names()`
- [packages/app-commands/src/commands/shared/run_paths.rs#L19](../../packages/app-commands/src/commands/shared/run_paths.rs#L19) and [packages/app-commands/src/commands/shared/run_paths.rs#L251-L258](../../packages/app-commands/src/commands/shared/run_paths.rs#L251-L258): `CDS_SELECTION_KEY` and `fn listed_cds_names()`

When the key or the shape of the setting changes, a change of one copy leaves the run paths and the amino-acid node data reading different CDS lists.

## Fix

Keep one reader of the CDS selection, in the module that owns the run settings, and call it from both places.

## Validation

`just test-rs -E 'package(app-commands)'` passes, and `rg 'CDS_SELECTION_KEY' packages` finds one definition.
