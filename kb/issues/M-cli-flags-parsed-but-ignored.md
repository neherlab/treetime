# CLI flags parsed but ignored

Several command-line flags are accepted by the parser, validated, and then never read. A user who sets one of them gets no error and no effect. This is the same class of defect as a silently wrong result: the command reports success while ignoring part of its input.

Two mechanisms hide the flags:

- **Resolved argument field never read**: the command's resolved argument struct (`Treetime<Command>Args` in `packages/app-cli/src/commands/<command>/args.rs`) stores the value, and no code reads it. The compiler detects these because the command modules are crate-private; each field carries `#[expect(dead_code, reason = "...")]` that points here, so the expectation fails the build's lint check as soon as a flag is wired
- **Core configuration field never read**: the command copies the value into a public core configuration struct, and the core never reads that field. The compiler cannot detect these because the core field is public

## Resolved argument fields never read

Flags marked *hidden* are accepted but not listed in `--help`.

| Command     | Flags                                                                                                                                                                                                              |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `ancestral` | `--zero-based`, `--report-ambiguous`, `--aa` (hidden), `--marginal` (hidden), `--custom-gtr` (hidden)                                                                                                              |
| `clock`     | `--alignment`/`--aln`, `--model`/`--gtr`, `--model-params`/`--gtr-params`, `--branch-length-mode`, `--method-anc`, `--prune-short`, `--seed`, `--clock-filter-method` (hidden), `--plot-rtt` (hidden), `--prune-outliers` (hidden) |
| `mugration` | `--seed`                                                                                                                                                                                                           |
| `timetree`  | `--tip-labels`, `--no-tip-labels`, `--n-iqd`, `--aa` (hidden), `--custom-gtr` (hidden), `--clock-filter-method` (hidden), `--greedy-resolve` (hidden), `--stochastic-resolve` (hidden)                            |

Tracked elsewhere, with their own `expect` reasons:

- `--vcf-reference` in `ancestral`, `clock`, and `timetree`: [M-io-vcf-input-output-unimplemented.md](M-io-vcf-input-output-unimplemented.md)
- `--method-anc` in `timetree`: [M-timetree-method-anc-ignored.md](M-timetree-method-anc-ignored.md)
- every `homoplasy` flag: [H-homoplasy-command-unimplemented.md](H-homoplasy-command-unimplemented.md)

## Core configuration fields never read

The `timetree` command copies these flags into `TimetreeParams` in `packages/treetime/src/timetree/pipeline.rs`, and the timetree pipeline never reads the fields:

- `--keep-polytomies` (`keep_polytomies`)
- `--report-ambiguous` (`report_ambiguous`)

## Impact

- `--prune-short` in `clock` leaves short branches in the tree
- `--seed` in `clock` and `mugration` does not make runs reproducible, because the value never reaches a random number generator
- `--model` and `--model-params` in `clock` do not change the substitution model
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
