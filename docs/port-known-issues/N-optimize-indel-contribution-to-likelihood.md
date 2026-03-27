# Branch length likelihood does not account for indels

Branch length optimization considers only nucleotide substitutions when computing the likelihood function. Insertions and deletions (indels) are ignored. A branch with zero substitutions but one or more indels is assigned zero length, even though the indel represents genuine evolutionary change.

## Current behavior

The per-edge likelihood is computed as a product over alignment positions:

$$L(t) = \prod_i \sum_{ab} s^i_a \, e^{Q t}_{ab} \, r^i_b$$

where $s^i$ and $r^i$ are the child-side and parent-side messages at position $i$, and $e^{Qt}$ is the substitution model transition matrix. The product runs over canonical (non-gap, non-ambiguous) positions only. Gap positions are excluded entirely via `edge_effective_length()`.

This means:

- `get_coefficients()` at [packages/treetime/src/commands/optimize/optimize_dense.rs#L59-L67](../../packages/treetime/src/commands/optimize/optimize_dense.rs#L59-L67) extracts eigenvalue-space coefficients from substitution messages only
- `evaluate_dense_contribution()` and `evaluate_sparse_contribution()` sum log-likelihoods over substitution-only site contributions
- `is_zero_branch_optimal()` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L192-L214](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L192-L214) tests the derivative at $t = 0$ using substitution contributions only
- `initial_guess_mixed()` at [packages/treetime/src/commands/optimize/optimize_unified.rs#L300-L328](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L300-L328) counts substitutions only via `edge_subs().len()`

The `DenseEdgePartition` struct at [packages/treetime/src/representation/payload/dense.rs#L36-L42](../../packages/treetime/src/representation/payload/dense.rs#L36-L42) stores `indels: Vec<InDel>`, and the prune command's `merge_sibling_pair()` preserves indels during topology changes, but no optimization code reads the indel field.

## Expected behavior

A branch carrying indels but no substitutions should have a positive optimal length proportional to the indel content. The likelihood function should include an indel contribution that penalizes zero-length branches when indels are present.

## Scientific background

Most phylogenetic software (RAxML, IQ-TREE, PhyML, BEAST) also ignores indels in the likelihood computation, treating gaps as missing data. TreeTime v0 follows this convention.

Incorporating indels into the likelihood requires modeling the indel process. Common approaches:

- **Affine gap penalty model**: add a fixed log-likelihood contribution per indel event (opening cost) plus a per-position extension cost. Simple to implement but not probabilistic.
- **TKF91/TKF92 models**: probabilistic indel models that extend the substitution likelihood with birth-death processes for insertions and deletions. Computationally expensive and rarely used in practice.
- **Empirical indel rate**: count indel events per branch, estimate an indel rate $\mu_{\text{indel}}$, and add $\log P(\text{n indels} \mid \mu_{\text{indel}} \cdot t)$ to the branch log-likelihood. A Poisson model $P(k \mid \lambda) = e^{-\lambda} \lambda^k / k!$ with $\lambda = \mu_{\text{indel}} \cdot t$ is the simplest option.

The empirical Poisson approach is the most practical for TreeTime's use case. The indel rate can be estimated from the alignment (total indel events / total branch length), and the per-branch contribution is a single scalar added to the log-likelihood sum.

## Proposed solution

Add an optional indel contribution to the per-edge log-likelihood:

1. Count indel events on each edge (already stored in `DenseEdgePartition.indels` and `SparseEdgePartition.indels`)
2. Estimate a global indel rate $\mu_{\text{indel}}$ from the tree (total indel events / total branch length / alignment length)
3. For each edge with $k$ indels and branch length $t$: add $\log P(k \mid \mu_{\text{indel}} \cdot L \cdot t)$ to the edge log-likelihood, where $L$ is the alignment length
4. Include the Poisson derivative contributions in the Newton step:
   - $\frac{d}{dt} \log P = (k / t) - \mu_{\text{indel}} \cdot L$
   - $\frac{d^2}{dt^2} \log P = -k / t^2$

This integrates into `OptimizationMetrics` by adding the indel terms to `log_lh`, `derivative`, and `second_derivative` after the substitution contributions.

For `is_zero_branch_optimal()`: if any edge has indels, the derivative at $t = 0$ from the indel Poisson term is $+\infty$ (since $k / t \to \infty$), which forces the derivative positive and prevents zero-length assignment. This is the desired behavior.

For `initial_guess_mixed()`: include indel events in the numerator (weighted by an indel-to-substitution rate ratio) or keep the substitution-only formula as an initial estimate and let Newton optimization find the correct length.

## Impact

Low for most datasets. Indels are rare compared to substitutions in typical viral phylogenetics (TreeTime's primary use case). The effect is most visible for:

- Trees with long branches spanning insertion/deletion events (bacterial phylogenetics, gene family evolution)
- Alignments with many gap characters relative to substitutions
- Branches where the only signal of divergence is an indel

## v0 handling

v0 ignores indels in the likelihood, same as current v1. This is a design-doc feature, not a v0 port gap.
