# Timetree sequence-evolution model text output lacks a replacement decision

TreeTime v0 writes `sequence_evolution_model.txt`; v1 writes `model.json` when a model is available. No approved decision establishes that the JSON output preserves the text file's information and workflow interoperability or intentionally replaces its public contract.

Compare both schemas and downstream uses. Implement the text contract unless an approved parity decision establishes `model.json` as its replacement.

Locations: v0 writes the text file in `packages/legacy/treetime/treetime/wrappers.py#L512` (with the tree suffix) and `packages/legacy/treetime/treetime/wrappers.py#L635`; v1 writes the model through `fn write_gtr_output()` [packages/app-commands/src/commands/shared/gtr_output.rs#L35](../../packages/app-commands/src/commands/shared/gtr_output.rs#L35). No file in `kb/decisions/` names `sequence_evolution_model`.

## Proposed implementation plan

1. Compare the content: list the fields of v0 `sequence_evolution_model.txt` (the `str()` of the GTR model, the substitution rate, and the equilibrium frequencies) against the fields of v1 `model.json`, and search the augur and nextstrain clones in `.repos/` for readers of the text file
2. When a field of the text file is missing from `model.json`, add it to the JSON output in `fn write_gtr_output()`
3. Without an approved replacement decision: add a writer of `sequence_evolution_model.txt` next to `fn write_gtr_output()` that renders the same model value in the v0 layout, and call it from the `timetree` output path. Test it in `packages/app-commands/src/commands/shared/__tests__/test_gtr_output.rs` with the text that v0 writes for the same model parameters as the oracle (an `indoc!` literal from a v0 run on `data/flu/h3n2/20`)
4. Validate with `just r treetime timetree` on `data/flu/h3n2/20` and a diff against the v0 file of the same run
