# Timetree with input branch lengths aborts on a point-by-point division

`treetime timetree --branch-length-mode=input` fails on every dataset the smoke matrix covers:

```
Error:
   0: Cannot divide point by point: operation not well-defined

Location: packages/treetime-distribution/src/distribution_ops/divide.rs:33
```

## Symptom and reproduction

```bash
treetime timetree --tree=data/flu/h3n2/20/tree.nwk --dates=data/flu/h3n2/20/metadata.tsv \
  --aln=data/flu/h3n2/20/aln.fasta.xz --branch-length-mode=input --output-all=<dir>
```

The same abort occurs on `data/ebola/20` and on `data/ebola/362` with `--keep-root` added. These are the three `--branch-length-mode=input` cases in `dev/run-smoke-tests` (`flu/h3n2/20/branch-input`, `ebola/20/branch-input`, `ebola/362/input-keeproot`), and all three fail.

## Impact and scope

- `--branch-length-mode=input` is unusable for `timetree`: the command produces no output at all. The alternative mode (`marginal`, the default) completes on the same datasets.
- The three cases are excluded from every before/after output comparison, so this region of the timetree pipeline carries no regression coverage.

## Mechanism

The time-marginal forward pass divides a node's composite distribution by the message it received, to form the cavity distribution. `fn distribution_division()` [packages/treetime-distribution/src/distribution_ops/divide.rs#L33](../../packages/treetime-distribution/src/distribution_ops/divide.rs#L33) rejects a `Point`/`Range` divisor outright: only division by a `Function` is defined.

With input branch lengths, the per-edge branch-length distribution is not built from a clock model, so a message can stay a degenerate `Point` (a single time with all the mass) instead of becoming a sampled `Function`. Dividing by that point has no well-defined quotient, and the pass aborts rather than treating the degenerate case.

Deciding the intended behavior is part of the fix: either the input-length path must produce non-degenerate messages, or division by a point divisor needs a defined result (for example, an empty quotient outside the point and an undefined-but-guarded quotient at it).

## Related issues

- [M-distribution-normalization-erases-errors.md](M-distribution-normalization-erases-errors.md)
- [M-timetree-internal-dates-missing-input-bl.md](M-timetree-internal-dates-missing-input-bl.md)
- [N-timetree-gtr-json-missing-input-bl.md](N-timetree-gtr-json-missing-input-bl.md)
