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

## Status

`fix(grid): follow the hard-approach law in a cell that ends at zero probability` on branch `claude-hard-cell-interp`
removes the cause traced below for a simulated reproducer (branch `claude-meanfield-prior`): `data/dengue/500` now
completes. `data/lassa/L/200` gets past this point already on `claude-work` and panics later converting a negative
date ([H-timetree-negative-numdate-panics-in-date-string-conversion.md](H-timetree-negative-numdate-panics-in-date-string-conversion.md)),
with or without the fix. The `dev/smoke.toml` declaration of the dengue case as an expected failure is not updated.

## Root cause

In a run whose rate had collapsed to $1.67 \cdot 10^{-5}$ (validation suite, `piecewise` rep09, $L = 10^3$,
$\mu = 3 \cdot 10^{-4}$, `--reroot least-squares --coalescent-opt --clock-filter 0`, mean-field coalescent prior), one
grid cell was 2.32 years. The parent's posterior stored negative log `+inf` only at its last grid point, the hard
boundary, and its peak in the point before. `GridFn::interpolate_at` interpolated negative logs linearly, so the whole
cell next to the boundary read as `+inf`, and a cell starting at `+inf` gave NaN. The cavity built on the same points
read the peak, a rounding error inside that cell, as `+inf`; resampling it for the convolution produced NaN;
`min_or` (ndarray-stats `min()` fails on NaN) made the operand look massless, and `distribution_convolution_fine`
returned `Empty`. The node's posterior was empty and its 45 descendants inherited an empty parent. A cell with an
infinite end now follows the hard-approach law of that boundary, anchored at its finite end.

## Mechanism

The division fails when one operand of a message division has empty support. An empty posterior arises when the message from the rest of the tree and a node's date constraint have disjoint support, the condition described in [H-timetree-backward-pass-nan-time-distribution.md](H-timetree-backward-pass-nan-time-distribution.md). The source of the empty operand has not been traced.

## Related issues

- [H-timetree-backward-pass-nan-time-distribution.md](H-timetree-backward-pass-nan-time-distribution.md)
- [H-timetree-input-branch-lengths-abort-on-point-division.md](H-timetree-input-branch-lengths-abort-on-point-division.md)
