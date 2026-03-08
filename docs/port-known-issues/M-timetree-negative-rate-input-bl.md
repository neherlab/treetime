# Negative clock rate silently stored in input-BL mode

`--branch-length-mode=input` bypasses the `debug_assert!(clock_rate >= 0.0)` in
[`branch_length_likelihood.rs#L38`](../../packages/treetime/src/commands/timetree/inference/branch_length_likelihood.rs#L38)
because the marginal solver path is skipped. The negative rate is stored in
`timetree.json` and used to compute dates, producing invalid output.

## Affected datasets

- rsv/a/20: rate = -0.916, 4/39 annotations
- dengue/500: rate = -0.224, 49/912 annotations

These are the same datasets that crash with
[Negative clock rate causes assertion failure](H-timetree-negative-clock-rate-assertion.md)
in default mode.

## Root cause

The clock regression still computes a negative rate for datasets with poor
temporal signal. In default mode, this hits the assertion in
`branch_length_likelihood.rs`. In `--branch-length-mode=input`, the code path
that would hit the assertion is skipped. The negative rate is stored in
`timetree.json` and the inverted time axis produces degenerate distributions
with very few valid annotations.

## Repro

```bash
./dev/docker/run ./dev/dev r treetime -- timetree --clock-filter=0 --branch-length-mode=input \
  --tree=data/rsv/a/20/tree.nwk --dates=data/rsv/a/20/metadata.tsv \
  --name-column=strain --date-column=num_date \
  --outdir=tmp/repro-neg-rate-bl data/rsv/a/20/aln.fasta.xz
jq -r '.clock_rate' tmp/repro-neg-rate-bl/timetree.json
# -0.9155515711664515
```

## Related issues

- [Negative clock rate causes assertion failure](H-timetree-negative-clock-rate-assertion.md)
  same root cause, different code path
