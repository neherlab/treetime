# Chapter 6: Zero-length branches -- detection and pruning

[Back to index](_index.md) | Previous: [Chapter 5: Branch length optimization](5-branch-length-optimization.md) | Next: [Chapter 7: Polytomy resolution](7-polytomy-resolution.md)

## What is a zero-length branch?

A **zero-length branch** is an edge in the phylogenetic tree with branch length equal to (or indistinguishable from) zero. This means the parent and child nodes are evolutionarily identical -- no substitutions occurred along that lineage, or substitutions that occurred exactly canceled out.

```
      root
     /    \
    X      C        Branch root-X has length 0.0001
   / \               (near-zero: ~1 substitution per 10,000 sites)
  A   B
```

When a branch length is near zero relative to the alignment length, two interpretations exist:

1. **Real short branch.** A and B genuinely diverged from C very recently. Node X existed briefly before splitting into A and B.
2. **Artifact of binary resolution.** The true tree is a three-way split. The tree builder created X because it only produces binary trees. The branch root-X has near-zero length because X should not exist.

Telling these apart is the central challenge. The detection criterion must be conservative enough to keep real short branches while collapsing artifacts.

## Why standard ML cannot produce exact zeros

Branch lengths are constrained to `[0, infinity)`. When the true branch length is zero, the ML estimate does not land exactly at zero -- it lands at a tiny positive value. This is the **boundary problem** in statistical estimation (Self and Liang 1987).

The intuition: even when parent and child have identical sequences, a tiny positive branch length `t = epsilon` fits the data slightly better than `t = 0`. The transition probability matrix at `t = epsilon` has small off-diagonal entries that accommodate any residual uncertainty in the marginal reconstruction. With more data, the estimate approaches zero but never reaches it.

This is the branch-length manifestation of the **star tree paradox** (Steel and Penny 2000): ML always prefers a resolved tree over the true star topology, because near-zero internal branches add tiny likelihood improvements at negligible parameter cost.

Self and Liang (1987) showed that when the true parameter lies on the boundary, the likelihood ratio test statistic follows a mixture of chi-squared distributions (`0.5 * chi^2(0) + 0.5 * chi^2(1)` for a one-sided boundary), not the standard chi-squared. This correction is used by the approximate likelihood ratio test (aLRT) for branch support (Anisimova and Gascuel 2006).

## Detection approaches

### Approach 1: bootstrap and aLRT

The **phylogenetic bootstrap** (Felsenstein 1985) resamples alignment columns, rebuilds trees, and reports the fraction containing each clade. Low bootstrap support (< 70%) indicates an unstable branch -- a candidate for zero length.

The **approximate likelihood ratio test** (Anisimova and Gascuel 2006) is faster: for each branch, compare the best tree against the best NNI rearrangement. The test statistic `2 * (log L_best - log L_second)` uses the Self-Liang mixture distribution. Anisimova et al. (2011) found aLRT and aBayes offer the highest statistical power while being orders of magnitude faster than bootstrap.

### Approach 2: information criteria

**AIC** (Akaike 1974): `-2 * log L + 2 * k`. **BIC** (Schwarz 1978): `-2 * log L + log(n) * k`. A polytomy has one fewer parameter than the resolved tree. If the likelihood gain from resolution does not justify the extra parameter, the criterion prefers the polytomy.

### Approach 3: regularization

The **adaptive LASSO** (Zhang, Dinh, and Matsen 2021; based on Tibshirani 1996 and Zou 2006) adds an L1 penalty `lambda * sum(w_i * |b_i|)` to the likelihood, driving truly zero branches to exactly zero. The adaptive weights `w_i = 1/|b_hat_i|^gamma` ensure branches with small initial estimates receive heavy penalties. This achieves **oracle properties**: the estimator correctly identifies zero branches with probability approaching 1.

## What TreeTime does

### v0: heuristic pruning inside the optimization loop

v0's `prune_short_branches()` uses a compound criterion:

```python
for node in self.tree.find_clades():
    if node.up is None or node.is_terminal():
        continue
    if (node.branch_length < 0.1 * self.one_mutation) and \
       (self.gtr.prob_t(node.up._cseq, node._cseq, 0.0,
                        pattern_multiplicity=self.data.multiplicity(mask=node.mask)) > 0.1):
        # collapse: reparent children to grandparent
        node.up.clades = [k for k in node.up.clades if k != node] + node.clades
        for clade in node.clades:
            clade.up = node.up
```

**Condition 1** (`branch_length < 0.1 / L`): the branch is shorter than 10% of one expected substitution. Fast filter.

**Condition 2** (`prob_t(parent, child, 0) > 0.1`): the GTR model evaluates the probability that parent and child sequences could have zero evolutionary distance. If the sequences differ at positions where the model assigns low probability to identity, the branch is kept despite being short.

The threshold 0.1 is heuristic. The 10% probability cutoff is conservative -- most identical sequence pairs at zero distance have probability near 1.0.

**Where it runs:** Inside each iteration of `optimize_tree()` (joint mode). Once after `optimize_tree_marginal()` (marginal mode). See [Chapter 9](9-iteration-loop.md).

v0 code: [`packages/legacy/treetime/treetime/treeanc.py#L1475-L1496`](../../../packages/legacy/treetime/treetime/treeanc.py#L1475-L1496)

### v1: derivative-sign detection without pruning

v1's `is_zero_branch_optimal()` evaluates the derivative of log-likelihood at `t = 0`:

```
d/dt log L(0) = sum_i [(sum_c k_{ic} * lambda_c) / (sum_c k_{ic})]
```

If this sum is negative, the likelihood decreases as t increases from zero, making zero a local maximum. The function returns true when:

1. Every site has positive, finite likelihood at t=0 (no degenerate coefficients)
2. The total derivative is negative and finite

This is the correct criterion for local optimality at the boundary. For models with unimodal per-edge likelihoods (JC69, F81), this is also the global optimum. For models with multiple local maxima (K2P, HKY, GTR), the shortcut certifies a local boundary maximum only. No heuristic thresholds needed for the detection itself.

**The gap:** v1 sets `branch_length = 0.0` when zero is optimal but does **not** collapse the edge. The node stays in the tree. Zero-length edges accumulate across iterations, wasting computation and preventing polytomy resolution. This is tracked as `M-optimize-no-topology-cleanup-in-loop` in the known issues ledger.

v1 code: [`packages/treetime/src/commands/optimize/optimize_unified.rs#L192-L214`](../../../packages/treetime/src/commands/optimize/optimize_unified.rs#L192-L214)

### Comparison

| Aspect                           | v0                                                    | v1                                                                        |
| -------------------------------- | ----------------------------------------------------- | ------------------------------------------------------------------------- |
| Detection criterion              | Branch length < threshold AND probability > threshold | Derivative of log-likelihood at zero < 0                                  |
| Thresholds                       | 0.1 \* one_mutation, 0.1 probability (heuristic)      | None (local optimality condition; exact for unimodal models such as JC69) |
| Action taken                     | Collapses edge, reparents children                    | Sets branch length to 0.0, keeps edge                                     |
| Creates polytomies               | Yes                                                   | No                                                                        |
| Requires reconstructed sequences | Yes (for prob_t)                                      | No (uses eigendecomposition coefficients)                                 |

## The composition problem when collapsing edges

When an edge is collapsed, mutations on the two merging edges must be **composed**, not unioned. The Markov semigroup property (Chapter 2) requires:

```
A->G on edge1, G->T on edge2  =>  net A->T  (composition)
A->G on edge1, G->A on edge2  =>  no change (cancellation)
```

Set-union would produce `{A->G, G->T}` (two entries at one position) or `{A->G, G->A}` (false double mutation). Both are incorrect.

v1's `collapse_sparse_edge()` in [`packages/treetime/src/commands/prune/run.rs`](../../../packages/treetime/src/commands/prune/run.rs) uses `iterator_union()` from `treetime_utils::iterator::union` -- the wrong operation. A correct composition operation does not exist yet and must be implemented before topology cleanup can be integrated into the optimization loop. Tracked as `M-prune-collapse-uses-union-not-composition`.

## References

- Self, S. G., and K.-Y. Liang. 1987. "Asymptotic Properties of Maximum Likelihood Estimators under Nonstandard Conditions." _JASA_ 82(398):605-610. https://doi.org/10.1080/01621459.1987.10478472
- Steel, M., and D. Penny. 2000. "Parsimony, Likelihood, and the Role of Models." _Mol. Biol. Evol._ 17(6):839-850. https://doi.org/10.1093/oxfordjournals.molbev.a026364
- Felsenstein, J. 1985. "Confidence Limits on Phylogenies." _Evolution_ 39(4):783-791. https://doi.org/10.1111/j.1558-5646.1985.tb00420.x
- Anisimova, M., and O. Gascuel. 2006. "Approximate Likelihood-Ratio Test for Branches." _Syst. Biol._ 55(4):539-552. https://doi.org/10.1080/10635150600755453
- Anisimova, M., et al. 2011. "Survey of Branch Support Methods." _Syst. Biol._ 60(5):685-699. https://doi.org/10.1093/sysbio/syr041
- Akaike, H. 1974. "A New Look at the Statistical Model Identification." _IEEE Trans. Autom. Control_ 19(6):716-723. https://doi.org/10.1109/TAC.1974.1100705
- Schwarz, G. 1978. "Estimating the Dimension of a Model." _Ann. Stat._ 6(2):461-464. https://doi.org/10.1214/aos/1176344136
- Tibshirani, R. 1996. "Regression Shrinkage and Selection Via the Lasso." _JRSS:B_ 58(1):267-288. https://doi.org/10.1111/j.2517-6161.1996.tb02080.x
- Zou, H. 2006. "The Adaptive Lasso and Its Oracle Properties." _JASA_ 101(476):1418-1429. https://doi.org/10.1198/016214506000000735
- Zhang, C., V. C. Dinh, and F. A. Matsen IV. 2021. "Nonbifurcating Phylogenetic Tree Inference via the Adaptive LASSO." _JASA_ 116(534):858-873. https://doi.org/10.1080/01621459.2020.1778481
- Lewis, P. O., M. T. Holder, and K. E. Holsinger. 2005. "Polytomies and Bayesian Phylogenetic Inference." _Syst. Biol._ 54(2):241-253. https://doi.org/10.1080/10635150590924208
