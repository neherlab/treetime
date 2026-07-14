# Match v0 terminal-edge reroot variance

Split the complete terminal-edge variance, including the leaf offset, between the parent and child sides as v0 does.

## Required changes

- Compute the full terminal variance once.
- Apply fractions $x$ and $1-x$ to the two sides of the candidate root split.
- Use the same convention in clock, timetree, and optimize reroot paths.
- Do not add a v1 divergence decision without explicit approval.

## Validation

- Side variances at $x=0$, $x=0.5$, and $x=1$.
- Two-tip and near-terminal optimum golden masters against v0.
- Shared clock/optimize result for the same objective and inputs.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/N-reroot-leaf-variance-offset-diverges-from-v0.md](../issues/N-reroot-leaf-variance-offset-diverges-from-v0.md)
- [kb/issues/N-reroot-missing-v0-golden-master.md](../issues/N-reroot-missing-v0-golden-master.md)
