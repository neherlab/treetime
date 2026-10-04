# CLI flags parsed but ignored

Several command-line flags are accepted by the parser, validated, and then never read. A user who sets one of them gets no error and no effect. This is the same class of defect as a silently wrong result: the command reports success while ignoring part of its input.

The command's resolved argument struct (`Treetime<Command>Args` in `packages/app-commands/src/commands/<command>/args.rs`) stores the value, and no code reads it. The compiler detects these fields because the command modules are crate-private; each field carries `#[expect(dead_code, reason = "...")]` that points here, so the expectation fails the build's lint check as soon as a flag is wired.

## Resolved argument fields never read

Flags marked _hidden_ are accepted but not listed in `--help`.

| Command     | Flags                                                                                                                                                                                                                                                        |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `ancestral` | `--model-params`/`--gtr-params`, `--zero-based`, `--aa` (hidden), `--marginal` (hidden), `--custom-gtr` (hidden)                                                                                                                                             |
| `clock`     | `--model`/`--gtr`, `--model-params`/`--gtr-params`, `--branch-length-mode`, `--method-anc`, `--prune-short`, `--clock-filter-method` (hidden), `--plot-rtt` (hidden), `--prune-outliers` (hidden)                           |
| `timetree`  | `--model-params`/`--gtr-params`, `--tip-labels`, `--no-tip-labels`, `--n-iqd`, `--keep-polytomies`, `--zero-based`, `--aa` (hidden), `--custom-gtr` (hidden), `--clock-filter-method` (hidden), `--greedy-resolve` (hidden), `--stochastic-resolve` (hidden) |

`optimize` also accepts `--model-params`/`--gtr-params` and never reads it. The flag lives in the shared `ModelArgs` struct in `packages/app-commands/src/commands/shared/model.rs`, and no command reads its `model_params` field. `ModelArgs` derives `Serialize`, which reads every field, so no `dead_code` expectation marks the field. For example, `ancestral --model k80 --model-params kappa=5` and `--model-params kappa=0.2` write byte-identical GTR files with the default K80 rates.

`clock` and `timetree` accept `--date-format` (config key `date_format`, default `%Y-%m-%d`) and never read it. The flag lives in the shared `DateColumnArgs` struct in `packages/app-commands/src/commands/shared/metadata.rs`, which also derives `Serialize`, so no `dead_code` expectation marks the field. The date reader in `packages/treetime-io/src/dates_csv.rs` parses every date cell with `DateParserOptions::default()`, and `struct DateParserOptions` has no format field, so a string date in a format outside the built-in list stays unreadable whatever `--date-format` says.

Tracked elsewhere, with their own `expect` reasons:

- `--alignment`/`--aln` in `clock`: [M-clock-alignment-ignored.md](M-clock-alignment-ignored.md)
- `--vcf-reference` in `ancestral`, `clock`, and `timetree`: [M-io-vcf-input-output-unimplemented.md](M-io-vcf-input-output-unimplemented.md)
- `--method-anc` in `timetree`: [M-timetree-method-anc-ignored.md](M-timetree-method-anc-ignored.md)
- every `homoplasy` flag: [H-homoplasy-command-unimplemented.md](H-homoplasy-command-unimplemented.md)

## Impact

- `--prune-short` in `clock` leaves short branches in the tree
- `--model` and `--model-params` in `clock` do not change the substitution model
- `--model-params` in `ancestral`, `optimize`, and `timetree` does not change the parameters of the selected model, so a named model always uses its default parameters
- `--date-format` in `clock` and `timetree` does not change how dates are parsed, although its help text says it controls the parsing of string dates
- `--greedy-resolve` and `--stochastic-resolve` in `timetree` do not select a polytomy resolution strategy; see [kb/proposals/timetree-stochastic-polytomy-resolution.md](../proposals/timetree-stochastic-polytomy-resolution.md)

## Potential solutions

- Implement the documented behavior of a flag, following v0, and add an end-to-end test that the flag changes the output
- Remove a flag, together with its help text and every reference, when v1 does not support the behavior

## Recommendation

Decide each flag on its own against v0 behavior and its owning feature issue. The flags differ in scientific meaning, input requirements, and parity constraints, so one aggregate implementation ticket would bundle unrelated decisions.

## Ticket readiness

No aggregate ticket is ready. Create one focused ticket per flag after its disposition is decided.

## Related issues

- [M-cli-help-text-defects.md](M-cli-help-text-defects.md): help text for several of these flags
- [M-clock-filter-residual-parity.md](M-clock-filter-residual-parity.md): clock filter behavior, related to `--clock-filter-method` and `--n-iqd`
