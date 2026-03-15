# Optimize command hardcodes GTR to JC69

The optimize command ignores the `--model` CLI flag and uses JC69 for all partitions. This is inconsistent with the ancestral command, which wires `--model` through `get_gtr_sparse()` / `get_gtr_dense()` and supports all named models plus inference from data.

## Problem

`run_optimize()` ([packages/treetime/src/commands/optimize/run.rs#L34-L145](../../packages/treetime/src/commands/optimize/run.rs#L34-L145)) destructures `model_name` from args at line 39 but never uses it. Both sparse and dense partitions receive hardcoded JC69:

- Sparse: [run.rs#L81](../../packages/treetime/src/commands/optimize/run.rs#L81) `jc69(JC69Params::default())?`
- Dense: [run.rs#L105](../../packages/treetime/src/commands/optimize/run.rs#L105) `jc69(JC69Params::default())?`

Four FIXME comments mark this as known ([run.rs#L65](../../packages/treetime/src/commands/optimize/run.rs#L65), [run.rs#L77-L78](../../packages/treetime/src/commands/optimize/run.rs#L77-L78)).

The `--model` flag accepts all `GtrModelName` variants (Infer, JC69, K80, F81, HKY85, T92, TN93, Jtt92) via [args.rs#L39](../../packages/treetime/src/commands/optimize/args.rs#L39). An end-user running `optimize --model=hky85` gets JC69 results with no warning.

## Impact

JC69 assumes equal base frequencies and equal substitution rates. For datasets with compositional bias (e.g. GC-rich genomes, amino acid data), JC69 produces systematically wrong branch lengths: overestimating branches between sequences with similar composition and underestimating branches between sequences with different composition.

## Consistency with other commands

| Command   | `--model` wired? | GTR inference? | Code path                |
| --------- | ---------------- | -------------- | ------------------------ |
| ancestral | Yes              | Yes            | `get_gtr_sparse/dense()` |
| timetree  | No (same bug)    | No             | hardcoded JC69           |
| optimize  | No (this issue)  | No             | hardcoded JC69           |
| mugration | N/A              | Partial        | mugration-specific GTR   |

The ancestral command ([packages/treetime/src/commands/ancestral/run.rs#L130](../../packages/treetime/src/commands/ancestral/run.rs#L130), [run.rs#L172](../../packages/treetime/src/commands/ancestral/run.rs#L172)) demonstrates the correct pattern: call `get_gtr_sparse(model_name, partition, &graph)` after initial Fitch compression, then replace the dummy GTR. The optimize command should follow the same pattern.

## Related issues

- [GTR model selection not implemented in timetree](M-timetree-gtr-selection.md): same bug in the timetree command
- [Dead CLI flags in timetree](N-timetree-dead-cli-flags.md): `--gtr-params` also dead in timetree

## Proposed solution

Follow the ancestral command's pattern:

1. Create partitions with dummy JC69 (existing code)
2. Run Fitch compression for sparse / `initialize_marginal()` for dense (existing code)
3. Call `get_gtr_sparse(model_name, partition, &graph)` or `get_gtr_dense()` to get the real GTR (replaces hardcoded JC69 with actual model dispatch)
4. Replace partition GTR (existing code, just with real model)
5. For `--model=infer` on dense: run a first marginal pass with dummy GTR to populate profiles, then infer GTR, then re-run (matching ancestral command's two-pass approach at [packages/treetime/src/commands/ancestral/run.rs#L168-L178](../../packages/treetime/src/commands/ancestral/run.rs#L168-L178))

The infrastructure (`get_gtr_by_name()`, `get_gtr_sparse()`, `get_gtr_dense()`, `infer_gtr_sparse()`, `infer_gtr_dense()`) already exists. The fix wires the existing `model_name` argument through to these functions.
