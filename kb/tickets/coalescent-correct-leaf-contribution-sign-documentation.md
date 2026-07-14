# Correct coalescent leaf contribution sign documentation

Document the existing $-I(t_{\mathrm{leaf}})$ value as a node-grouped negative-log term and add a direct branch-sum test that prevents an incorrect sign change.

## Required changes

1. Rewrite the leaf contribution documentation in `packages/treetime/src/coalescent/contributions.rs` so it distinguishes an isolated branch-survival cost from the node-grouped whole-tree decomposition.
2. State that `DistributionNegLog` stores negative-log terms even though an individual grouped term can be negative; remove the claim that the type contains log probabilities.
3. Preserve the numerical leaf formula $-I(t_{\mathrm{leaf}})$.
4. Add a direct unit test for a fixed tree. If $c$ and $p$ denote a branch's child and parent, respectively, assert that the grouped survival terms equal

   $$
   \sum_{(c,p)} \left[I(t_p)-I(t_c)\right].
   $$

5. Cite the equivalent v0 leaf and root grouping in the test oracle comment.

## Validation

- Run the focused coalescent contribution tests.
- Run the full lint and test suite.

## Related issues

- Source: [kb/issues/N-coalescent-leaf-contribution-comments-misstate-sign.md](../issues/N-coalescent-leaf-contribution-comments-misstate-sign.md)
- Related: [kb/issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md](../issues/M-timetree-coalescent-missing-leaf-and-root-contributions.md)
