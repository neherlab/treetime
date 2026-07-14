# Coalescent backward pass misses leaf and root contributions

V1 computes node-grouped leaf and internal coalescent contributions, but the backward pass omits leaf contributions and no root correction is constructed. This leaves the telescoped whole-tree coalescent prior incomplete [packages/treetime/src/timetree/inference/backward_pass.rs#L49-L53](../../packages/treetime/src/timetree/inference/backward_pass.rs#L49-L53) [packages/treetime/src/coalescent/contributions.rs#L31-L76](../../packages/treetime/src/coalescent/contributions.rs#L31-L76).

## Mathematical contract

Let $I(t)$ be the cumulative per-branch merger rate at time $t$, $m_j$ be the merger multiplicity of internal node $j$, and $\lambda(t)$ be the instantaneous total merger rate. The node-grouped negative-log terms are:

- internal non-root node: $m_j[I(t_j)-\log\lambda(t_j)]$;
- leaf: $-I(t_j)$;
- root correction: $+I(t_{\mathrm{root}})$ once, after combining the root's child messages.

For a branch from child time $t_c$ to parent time $t_p$, survival contributes $I(t_p)-I(t_c)$ to negative log-likelihood. Summing branches gives each non-root internal node $m_j I(t_j)$, each leaf $-I(t_j)$, and an additional $+I(t_{\mathrm{root}})$ because the root has no parent branch. V0 applies the same leaf and root terms [packages/legacy/treetime/treetime/clock_tree.py#L499-L530](../../packages/legacy/treetime/treetime/clock_tree.py#L499-L530).

## Impact

- The coalescent prior propagated through the tree is incomplete.
- Root and internal node-time posteriors can differ from v0.
- Constant-$T_c$, optimized-$T_c$, and skyline inference all consume the affected backward pass.

## Recommendation

Apply the existing $-I(t_{\mathrm{leaf}})$ contribution before each leaf message is convolved toward its parent. After combining all root child messages, add $+I(t_{\mathrm{root}})$ exactly once, evaluating the TBP-defined integral through the established calendar-to-TBP conversion. Validate the complete likelihood against a direct branch-sum calculation and v0.

## Related issues

- [H-timetree-coalescent-events-incomplete-after-topology-change.md](H-timetree-coalescent-events-incomplete-after-topology-change.md)
- [M-timetree-positional-likelihood-metric.md](M-timetree-positional-likelihood-metric.md)
