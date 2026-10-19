# The coalescent prior counts a leaf without a date only at its parent

A leaf without a date is a bad branch [packages/treetime/src/timetree/inference/bad_branches.rs#L7-L22](../../packages/treetime/src/timetree/inference/bad_branches.rs#L7-L22). Under a coalescent prior, time inference treats its lineage in three different ways:

- **Its own posterior has no coalescent factor.** The backward pass builds no subtree distribution for it: with no child messages and no date constraint the distribution is `Empty`, and `fn apply_leaf_coalescent_factor()` multiplies `Empty` by the leaf factor, which leaves it `Empty` [packages/treetime/src/timetree/inference/backward_pass.rs#L53-L59](../../packages/treetime/src/timetree/inference/backward_pass.rs#L53-L59), [packages/treetime/src/timetree/inference/backward_pass.rs#L135-L158](../../packages/treetime/src/timetree/inference/backward_pass.rs#L135-L158), [packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L15-L17](../../packages/treetime-distribution/src/distribution_ops/multiply_by_fn.rs#L15-L17). The forward pass then takes the forward message alone as its posterior [packages/treetime/src/timetree/inference/forward_pass.rs#L145-L150](../../packages/treetime/src/timetree/inference/forward_pass.rs#L145-L150). The leaf factor `node_contribution(key, 0, t)`, the survival term of the leaf's branch at its lower end, is missing from the posterior of the leaf date; a leaf has no merger term
- **The expected lineage counts leave it out.** `fn LineageExpectation::from_posterior()` gives a bad branch no survival function [packages/treetime/src/coalescent/mean_field.rs#L176-L184](../../packages/treetime/src/coalescent/mean_field.rs#L176-L184)
- **Its parent counts it.** The parent's factor uses the number of all its children, bad branches included [packages/treetime/src/timetree/inference/backward_pass.rs#L56](../../packages/treetime/src/timetree/inference/backward_pass.rs#L56), so the merger term and the survival terms of the parent count the lineage of the undated leaf [packages/treetime/src/coalescent/mean_field.rs#L94-L113](../../packages/treetime/src/coalescent/mean_field.rs#L94-L113)

[kb/decisions/timetree-uncertain-leaf-dates-are-inferred.md](../decisions/timetree-uncertain-leaf-dates-are-inferred.md) states that an undated leaf "is excluded from the coalescent". The parent's count contradicts that statement; the missing leaf factor agrees with it.

## v0

v0 behaves the same way: a leaf without date constraint and without child messages gets no subtree distribution and no merger contribution [packages/legacy/treetime/treetime/clock_tree.py#L711-L722](../../packages/legacy/treetime/treetime/clock_tree.py#L711-L722), and its posterior is the message from its parent [packages/legacy/treetime/treetime/clock_tree.py#L889-L891](../../packages/legacy/treetime/treetime/clock_tree.py#L889-L891). The lineage counts skip bad branches [packages/legacy/treetime/treetime/merger_models.py#L102-L105](../../packages/legacy/treetime/treetime/merger_models.py#L102-L105), while the parent's contribution uses `len(node.clades)` [packages/legacy/treetime/treetime/merger_models.py#L239-L250](../../packages/legacy/treetime/treetime/merger_models.py#L239-L250). The comment at the call says "set multiplicity of node to number of good child branches", but the call passes no multiplicity [packages/legacy/treetime/treetime/clock_tree.py#L735-L743](../../packages/legacy/treetime/treetime/clock_tree.py#L735-L743).

## Model consequence

> [!IMPORTANT]
> **Investigation required.** The size of the effect. The parent's merger term counts one more merging lineage than the lineage counts hold, and the date of the undated leaf ignores the survival of its branch. Measure on a dataset with a coalescent prior (`--coalescent-opt`) and a few dates removed: the dates and 90% intervals of the undated leaves and their parents, and the fitted $T_c$, with the parent counting only children that are not bad branches, against the current output.

## Fix direction

> [!IMPORTANT]
> **Decision required.** Treat an undated leaf in one way. Options: (a) exclude it fully: parents count only children that are not bad branches, which matches the stated decision and the comment in v0 but not the code of v0; (b) include it fully: the leaf factor enters the leaf's posterior in the forward pass and the leaf's lineage enters the expected lineage counts, which changes the stated decision. Both diverge from v0 output.
