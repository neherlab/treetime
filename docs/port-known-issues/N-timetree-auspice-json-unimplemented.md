# Timetree does not write auspice JSON output

The timetree command does not produce `auspice_tree.json`. v0 writes this file for visualization in Nextstrain's Auspice tool. The design document (`docs/algorithms/timetree-optional-bits.md:28`) ties auspice output to confidence intervals: "In the auspice json uncertainty can be specified."

## Current state

Auspice types and serialization exist in `packages/treetime-io/src/auspice_types.rs` and `packages/treetime-io/src/auspice.rs`. A convert command in `packages/treetime-cli/src/convert/auspice.rs` can produce auspice JSON from other tree formats. But the timetree command has no code path to write auspice output with inferred dates and confidence intervals.

The feature inventory lists `[ ] Auspice JSON (v0 writes auspice_tree.json)` under timetree output.

## v0 implementation

`packages/legacy/treetime/treetime/CLI_io.py` writes `auspice_tree.json` with node dates and confidence intervals embedded in the JSON structure. Auspice displays these as date range bars on the tree visualization.

## Impact

Nextstrain users who run timetree and want to visualize results in Auspice must manually convert the output. The confidence interval data (now implemented in v1) cannot reach Auspice without this output format.
