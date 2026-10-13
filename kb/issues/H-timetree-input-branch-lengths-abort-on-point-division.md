# Timetree with input branch lengths aborts on a point-by-point division

> [!IMPORTANT]
> **Investigation required.** The code path that turns an input branch length into a `Point` message is not traced. Find where v1 builds the branch-length distribution in input mode, and confirm that the abort disappears once it is built as v0 builds it (Fix direction below).

`treetime timetree --branch-length-mode=input` fails on every dataset the smoke matrix covers:

```
Error:
   0: Cannot divide point by point: operation not well-defined

Location: packages/treetime-distribution/src/distribution_ops/divide.rs:34
```

## Symptom and reproduction

```bash
treetime timetree --tree=data/flu/h3n2/20/tree.nwk --dates=data/flu/h3n2/20/metadata.tsv \
  --aln=data/flu/h3n2/20/aln.fasta.xz --branch-length-mode=input --output-all=<dir> \
  --output-selection=nwk,nexus,auspice,augur-node-data,clock-model
```

Input mode reconstructs no ancestral sequences, so `--output-all` skips the reconstructed FASTA, and `--output-reconstructed-nuc-fasta` fails before inference.

The same abort occurs on `data/ebola/20` and on `data/ebola/362` with `--keep-root` added. These are the `--branch-length-mode=input` cases of `dev/smoke` (the `branch-input`, `branch-input-all-outputs` and `input-keeproot` variants of `timetree`), and all of them fail.

The same abort occurs on the simulated data of [validation/README.md](../../validation/README.md).

## Impact and scope

- `--branch-length-mode=input` is unusable for `timetree`: the command produces no output at all. The alternative mode (`marginal`, the default) completes on the same datasets.
- These cases are excluded from every before/after output comparison, so this region of the timetree pipeline carries no regression coverage.

## Mechanism

The time-marginal forward pass divides a node's composite distribution by the message it received, to form the cavity distribution. `fn distribution_division()` [packages/treetime-distribution/src/distribution_ops/divide.rs#L31-L34](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L31-L34) rejects a `Point`/`Range` divisor outright: only division by a `Function` is defined.

With input branch lengths, the per-edge branch-length distribution is not built from a clock model, so a message can stay a degenerate `Point` (a single time with all the mass) instead of becoming a sampled `Function`. Dividing by that point has no well-defined quotient, and the pass aborts rather than treating the degenerate case.

## Fix direction

Match v0, which builds a non-degenerate branch-length distribution from the input length [packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102](../../packages/legacy/treetime/treetime/branch_len_interpolator.py#L64-L102):

- branch length below 0.05: Poisson likelihood of the number of mutations
- longer branches: Gaussian, with a variance corrected for saturation

## Validation

The `--branch-length-mode=input` cases of `dev/smoke` and of the validation suite complete, and the per-edge distribution in input mode matches v0's on a small tree.

## Related issues

- [M-timetree-internal-dates-missing-input-bl.md](M-timetree-internal-dates-missing-input-bl.md)
- [N-timetree-gtr-json-missing-input-bl.md](N-timetree-gtr-json-missing-input-bl.md)
