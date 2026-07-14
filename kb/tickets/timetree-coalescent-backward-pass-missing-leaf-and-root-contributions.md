# Apply complete coalescent leaf and root contributions

Apply the existing node-grouped leaf term and the exact root correction in timetree belief propagation.

## Required changes

1. Preserve the leaf helper's node-grouped negative-log value $-I(t_{\mathrm{leaf}})$, where $I(t)$ is cumulative per-branch merger rate.
2. Apply each leaf contribution before convolving the leaf message toward its parent.
3. After combining the root's child messages, add $+I(t_{\mathrm{root}})$ exactly once.
4. Evaluate $I$ in TBP coordinates using the established calendar-to-TBP conversion at both leaf and root boundaries.
5. Include all contributions in reported coalescent likelihoods and convergence state.

## Validation

- Assert the node-grouped sum equals the direct branch survival sum $\sum_{(c,p)}[I(t_p)-I(t_c)]$ plus internal merger-rate terms.
- Direct analytical trees with precise and uncertain leaf dates.
- Whole-tree v0 golden masters after explicit sign-space conversion.
- Root, leaf, and internal contribution sum against an independent direct likelihood.
- Full lint and test suite.

## Related issues

- Source: [kb/issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md](../issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md)
- Related: [kb/issues/H-timetree-coalescent-events-incomplete-after-topology-change.md](../issues/H-timetree-coalescent-events-incomplete-after-topology-change.md)
