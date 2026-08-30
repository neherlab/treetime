# timetree tip-state flags are accepted but not wired to output

## Symptom

`treetime timetree` accepts `--include-leaves`, `--impute-missing-data`, and the `--reconstruct-tip-states` alias, and threads them into `TimetreeParams` (`include_leaves`, `impute_missing_data`). No code in the timetree pipeline reads these fields, so the flags do not change any output. The predecessor `reconstruct_tip_states` field had the same defect: it was set from CLI args but never read.

## Impact and scope

A user who passes `--reconstruct-tip-states` (or the split flags) to `timetree` expecting imputed or emitted tip sequences gets no effect. The flags are inert. Scope is the `timetree` command only; the same flags on `ancestral` are fully wired.

## Root cause

`timetree` does not run the tip-emitting reconstruction path. Its node-data output derives tip sequences through its own serialization (`packages/treetime/src/commands/timetree/output/augur_node_data.rs`), which does not consult the tip-state parameters. The parameters exist on `TimetreeParams` (`packages/treetime/src/timetree/pipeline.rs`) but have no reader.

## Fix approach

Wire the timetree tip output through the same flag-aware reconstruction the `ancestral` command uses (`fn ancestral_reconstruction_marginal()` writes each partition's stored sequence, which the node-data serializer reads back), passing `include_leaves` and `impute_missing_data` from `TimetreeParams`. Then update the `--include-leaves` / `--impute-missing-data` / `--reconstruct-tip-states` row in `kb/features/timetree.md`.
