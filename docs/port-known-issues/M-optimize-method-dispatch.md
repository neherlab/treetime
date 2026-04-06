# Per-edge optimization: dispatch all 6 methods

## Problem

After the scaffolding, Brent, Newton-log, and step-clamping issues are implemented, all 6 inner-loop functions exist but only 5 are wired into `run_optimize_mixed`. The `NewtonLog` arm in the dispatch match still uses `todo!()`.

`BrentSqrt` and `BrentLog` were wired as part of M-optimize-method-brent-parameterized.

## Context

The 6 inner-loop functions and their locations (after M-optimize-method-scaffolding split the file):

| Variant      | Function            | Module        | Status |
| :----------- | :------------------ | :------------ | :----- |
| `Brent`      | `brent_inner`       | Brent module  | Wired  |
| `BrentSqrt`  | `brent_sqrt_inner`  | Brent module  | Wired  |
| `BrentLog`   | `brent_log_inner`   | Brent module  | Wired  |
| `Newton`     | `newton_inner`      | Newton module | Wired  |
| `NewtonSqrt` | `newton_sqrt_inner` | Newton module | Wired  |
| `NewtonLog`  | `newton_log_inner`  | Newton module | todo!  |

All share the same prelude in `run_optimize_mixed` (contribution collection, indel setup, zero-branch check, min_branch_length computation). The only difference: Newton variants need `evaluate_with_indels` to get metrics (derivatives), Brent variants need only the branch length and contributions.

## Approach

Replace `todo!()` arm for `NewtonLog` with a call to the new function. Add necessary import from the Newton module. The dispatch block is at `optimize_unified.rs` (currently `run_optimize_mixed`).

## Verification

`./dev/docker/run ./dev/dev t` -- existing tests pass (they use `Newton`, `NewtonSqrt`, `Brent`).

Manual smoke test for each new variant:

```bash
for method in brent-sqrt brent-log newton-log; do
  ./dev/docker/run ./dev/dev r treetime -- optimize \
    --opt-method=$method \
    --tree=data/flu/h3n2/20/tree.nwk \
    --aln=data/flu/h3n2/20/aln.fasta.xz \
    --outdir=tmp/opt-$method
done
```

## Dependencies

- Depends on: ~~M-optimize-method-brent-parameterized~~ (done), [M-optimize-method-newton-log](M-optimize-method-newton-log.md), [M-optimize-method-step-clamping](M-optimize-method-step-clamping.md)
- Depended on by: [M-optimize-method-tests](M-optimize-method-tests.md), [L-optimize-method-docs](L-optimize-method-docs.md)

## Cross-references

- Current dispatch: `packages/treetime/src/commands/optimize/optimize_unified.rs` (`run_optimize_mixed`)
