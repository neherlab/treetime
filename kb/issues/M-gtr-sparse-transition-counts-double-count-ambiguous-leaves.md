# Sparse transition counts probably count ambiguous leaf positions twice

`fn accumulate_sparse_transitions()` [packages/treetime/src/partition/marginal/sparse/count.rs#L66-L103](../../packages/treetime/src/partition/marginal/sparse/count.rs#L66-L103) sums the expected transitions of an edge over the explicit positions of `msg_to_parent` and over its fixed states, where each fixed state has the weight of its fixed `count`.

On a leaf edge, `msg_to_parent` is the leaf message (`fn leaf_message()` [packages/treetime/src/partition/marginal/sparse/backward.rs#L96-L119](../../packages/treetime/src/partition/marginal/sparse/backward.rs#L96-L119)). Its explicit positions are the ambiguous positions of the leaf. Its fixed profiles cover every determined state, ambiguous states included, each with the composition count of that state in the leaf sequence. An `R` position of a leaf therefore enters twice: once as an explicit position with the profile of `R`, and once in the fixed group of `R`, whose count includes it. Internal messages do not have this overlap, because their fixed counts exclude the positions that stay explicit.

> [!IMPORTANT]
> **Investigation required.** Confirm the double count with a test: a tree with one ambiguous leaf position, with the transition counts `nij` and the times per state `Ti` compared with a hand computation of the expected counts.

## Related issues

- [M-gtr-sparse-transition-counts-diverge-from-dense.md](M-gtr-sparse-transition-counts-diverge-from-dense.md): the sparse transition counts differ from the dense counts
