# Timetree aborts on dengue/500 and lassa/L/200 dividing by an empty distribution

`treetime timetree` with default options fails on `data/dengue/500` in round 1 and on `data/lassa/L/200` in round 2:

```
Error:
   0: When running round 1
   1: Timetree inference failed
   2: Cannot divide by empty distribution
```

The `rust` branch fails on both datasets with the same error, so the defect predates the current development branch.

## Symptom and reproduction

```bash
treetime timetree --tree=data/dengue/500/tree.nwk --dates=data/dengue/500/metadata.tsv \
  --aln=data/dengue/500/aln.fasta.xz --name-column=genbank_accession --output-all=<dir>
treetime timetree --tree=data/lassa/L/200/tree.nwk --dates=data/lassa/L/200/metadata.tsv \
  --aln=data/lassa/L/200/aln.fasta.xz --name-column=accession --output-all=<dir>
```

These are the `timetree/dengue/500/basic` and `timetree/lassa/L/200/basic` cases of `dev/smoke` (full tier).

## Impact and scope

- The command produces no output for these datasets.
- Both cases are declared expected failures in `dev/smoke.toml`, so a fix shows up as an unexpected pass until the declaration is removed.

## Mechanism

The division fails when one operand of a message division has empty support. An empty posterior arises when the message from the rest of the tree and a node's date constraint have disjoint support, the condition described in [H-timetree-backward-pass-nan-time-distribution.md](H-timetree-backward-pass-nan-time-distribution.md). The source of the empty operand has not been traced.

## Simulated reproducer without an abort

Validation suite, `piecewise` rep09, $L = 10^3$, $\mu = 3 \cdot 10^{-4}$ (no usable temporal signal), `--reroot least-squares
--coalescent-opt --clock-filter 0`: the root-to-tip regression after rerooting gives a rate of $1.67 \cdot 10^{-5}$
(true $3 \cdot 10^{-4}$). With the mean-field coalescent prior (branch `claude-meanfield-prior`) the marginal-likelihood
rate steps continue to $1.03 \cdot 10^{-5}$, about 100 years per substitution on a 1000-site alignment, and the forward
pass leaves 46 internal nodes with an empty posterior ("has an empty time distribution; no date was assigned"); the run
exits 0 and writes those nodes without a date. No node is flagged as a bad branch. With the frozen lineage counts the
rate stays at $1.65 \cdot 10^{-5}$ and every node gets a date.

## Root cause in the simulated reproducer

Traced with temporary diagnostics in the reproducer above; the first of the 46 undated nodes is the only one with a
non-empty parent posterior, the other 45 inherit an empty parent:

1. The parent's posterior has zero probability (negative log `+inf`) at its last two grid points, at the hard upper
   bound from its children; the cavity, the posterior divided by the node's upward message, keeps them.
2. `convolve_across_edge` calls `distribution_convolution_fine`, which resamples both operands to a common `dx`
   (`resample_dx`, here 2.3236 against the branch's 2.4635). `GridFn::interpolate_at`
   ([`grid_fn.rs#L237-L250`](../../packages/treetime-grid/src/grid_fn.rs#L237-L250)) computes `y0 + t * (y1 - y0)` on
   negative log values, which is NaN whenever `y0 = +inf`: 3 of the 416 resampled cavity values are NaN.
3. `to_peak_normalized_plain` takes the peak with `min_or`, whose ndarray-stats `min()` fails on any NaN; the default
   `+inf` makes the operand look massless, and the convolution returns `Empty` through `guarded_empty_result`.
4. The node's posterior is `Empty`, so its children have no parent posterior: all 45 are left undated in turn.

Fix direction: interpolation in negative log space treats a cell with an infinite end as having zero probability
(`+inf`) instead of producing NaN; whether `min_or` should also ignore NaN is a separate choice. Whether the aborts on
`data/dengue/500` and `data/lassa/L/200` share this cause is not verified.

## Related issues

- [H-timetree-backward-pass-nan-time-distribution.md](H-timetree-backward-pass-nan-time-distribution.md)
- [H-timetree-input-branch-lengths-abort-on-point-division.md](H-timetree-input-branch-lengths-abort-on-point-division.md)
