# CLI flags parsed but ignored

> [!IMPORTANT]
> **Decision required.** Each flag needs its own disposition against v0 behavior and its owning feature issue: implement the documented behavior, or remove the flag. See "Potential solutions".

Several command-line flags are accepted by the parser, validated, and then never read. A user who sets one of them gets no error and no effect. This is the same class of defect as a silently wrong result: the command reports success while ignoring part of its input.

The command's raw argument struct (`Treetime<Command>ArgsRaw` in `packages/app-commands/src/commands/<command>/args.rs`) stores the value. The command resolves it into `Treetime<Command>Args`, which keeps the raw struct and gives access to it through `Deref` (for example [packages/app-commands/src/commands/timetree/args.rs#L40-L84](../../packages/app-commands/src/commands/timetree/args.rs#L40-L84)). Neither the resolution nor the run reads these fields. The raw struct derives the deser `Serialize` and the clap parser, which read every field, so the compiler does not report these fields; the table below is the list to keep current when a flag is wired or removed.

## Argument fields never read

Flags marked _hidden_ are accepted but not listed in `--help`.

| Command     | Flags                                                                                                                                                                                                                                                        |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `ancestral` | `--model-params`/`--gtr-params`, `--zero-based`, `--aa` (hidden), `--marginal` (hidden), `--custom-gtr` (hidden)                                                                                                                                             |
| `clock`     | `--model`/`--gtr`, `--model-params`/`--gtr-params`, `--branch-length-mode`, `--method-anc`, `--prune-short`, `--clock-filter-method` (hidden), `--plot-rtt` (hidden), `--prune-outliers` (hidden)                                                            |
| `timetree`  | `--model-params`/`--gtr-params`, `--tip-labels`, `--no-tip-labels`, `--n-iqd`, `--keep-polytomies`, `--zero-based`, `--aa` (hidden), `--custom-gtr` (hidden), `--clock-filter-method` (hidden) |

`optimize` and `homoplasy` also accept `--model-params`/`--gtr-params` and never read it. The flag lives in the shared `ModelArgs` struct in `packages/app-commands/src/commands/shared/model.rs`, and no command reads its `model_params` field. The `--help` text and `docs/docs/reference.md` even give a usage example (`--model k80 --model-params kappa=0.2 pis=0.25,0.25,0.25,0.25`). Yet `ancestral --model k80 --model-params kappa=5` and `--model-params kappa=0.2` write byte-identical GTR files with the default K80 rates.

Tracked elsewhere, with their own `expect` reasons:

- `--alignment`/`--aln` in `clock`: [M-clock-alignment-ignored.md](M-clock-alignment-ignored.md)
- `--vcf-reference` in `ancestral`, `clock`, and `timetree`: [M-io-vcf-input-output-unimplemented.md](M-io-vcf-input-output-unimplemented.md)
- `--method-anc` in `timetree`: [M-timetree-method-anc-ignored.md](M-timetree-method-anc-ignored.md)

## Impact

- `--prune-short` in `clock` leaves short branches in the tree
- `--model` and `--model-params` in `clock` do not change the substitution model
- `--model-params` in `ancestral`, `homoplasy`, `optimize`, and `timetree` does not change the parameters of the selected model, so a named model always uses its default parameters

## Potential solutions

- Implement the documented behavior of a flag, following v0, and add an end-to-end test that the flag changes the output
- Remove a flag, together with its help text and every reference, when v1 does not support the behavior

## Recommendation

Decide each flag on its own against v0 behavior and its owning feature issue. The flags differ in scientific meaning, input requirements, and parity constraints, so one aggregate change would bundle unrelated decisions.

## Guard against new unread flags

No test checks that a parsed argument field is used. Add one test that fails when a field of a command argument struct is neither read by the command nor rejected during validation, so that a new flag cannot be added without an effect.

## Readiness

Implement each flag as a separate change after its disposition is decided. The guard test can land first, with the current unread fields as its known exceptions.

## Related issues

- [M-cli-help-text-defects.md](M-cli-help-text-defects.md): help text for several of these flags
- [M-clock-filter-residual-parity.md](M-clock-filter-residual-parity.md): clock filter behavior, related to `--clock-filter-method` and `--n-iqd`
