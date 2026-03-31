# Chapter 2: Gap treatment and counting models

[Back to index](_index.md) | Previous: [Chapter 1: Introduction](1-introduction.md) | Next: [Chapter 3: Single-residue birth-death models](3-single-residue.md)

## Gaps as missing data

Most phylogenetic ML software treats gap characters as missing data. At a gapped position, the partial likelihood vector is set to all 1.0 (uniform over states), contributing zero information to the site likelihood. The site likelihood at a gap position marginalizes over all possible states, which is equivalent to excluding that position from the substitution likelihood product.

### Why this is problematic

<a id="cite-1"></a>[Warnow 2012](https://doi.org/10.1371/currents.rrn1308) [[1](#ref-1)] constructed a four-taxon counterexample where ML phylogeny estimation treating indels as missing data is statistically inconsistent: the ML tree converges to the wrong topology as sequence length increases, because the missing-data treatment systematically underestimates distances between sequences with different gap patterns.

<a id="cite-27"></a>[Truszkowski and Goldman 2015](https://doi.org/10.1093/sysbio/syv089) [[27](#ref-27)] subsequently proved that ML IS consistent on gapped MSAs under broader (and more realistic) conditions: substitution rates $> 0$ on all edges and each site category observed infinitely often as the alignment grows. The Warnow counterexample relies on zero substitution rates on specific edges, which violates these conditions.

The consistency debate concerns tree topology estimation, not branch length optimization. Regardless of consistency, the practical consequence for branch length optimization remains: a branch with zero substitutions but one or more indels is assigned zero length, collapsing topology that the indel evidence supports. This is the direct motivation for TreeTime v1's Poisson indel contribution.

### Software using this approach

RAxML-NG, IQ-TREE, PhyML, BEAST (default), MrBayes, HyPhy, FastTree, TreeTime v0. In IQ-TREE, `STATE_UNKNOWN` sets the partial likelihood to `1.0` for all states (`alignment/pattern.cpp`). In libpll-2 (RAxML-NG), gap `-` maps to bitmask `1111` (all nucleotide states) in `src/maps.c`. In TreeTime v0, `GTR.state_pair()` uses `ignore_gaps=True` (default) to exclude gap positions from the compressed state-pair representation used by `optimal_t_compressed()`.

---

## Poisson indel count model (implemented in v1)

Models the number of indel events on each edge as an independent Poisson process with a single global rate. Each indel event contributes equally regardless of length or direction (insertion vs deletion). This is a TreeTime v1 design-doc feature with no counterpart in v0 or other phylogenetic software.

### Mathematical formulation

For a branch with $k$ indel events and length $t$, and global indel rate $\mu_{\text{indel}}$:

$$\ell_{\text{indel}}(t) = k \ln(\mu_{\text{indel}} t) - \mu_{\text{indel}} t - \ln(k!)$$

First and second derivatives with respect to $t$:

$$\frac{d\ell}{dt} = \frac{k}{t} - \mu_{\text{indel}}, \qquad \frac{d^2\ell}{dt^2} = -\frac{k}{t^2}$$

The global rate is estimated by the Poisson MLE over all edges:

$$\hat{\mu}_{\text{indel}} = \frac{\sum_e k_e}{\sum_e t_e}$$

where $k_e$ is the indel count on edge $e$ and $t_e$ is the current branch length estimate. The rate is re-estimated at each optimization round as branch lengths converge.

### Properties relevant to branch length optimization

The Poisson log-likelihood $\ell(t)$ is strictly concave in $t$ for $k > 0$ (the second derivative $-k/t^2 < 0$), guaranteeing a unique optimum at $t^* = k / \mu_{\text{indel}}$. At $t = 0$ with $k > 0$, $\ell \to -\infty$ and $d\ell/dt \to +\infty$, which forces the Newton optimizer away from zero. The indel contribution is additive with the substitution log-likelihood: the combined Newton step uses $\ell_{\text{sub}}(t) + \ell_{\text{indel}}(t)$ with summed derivatives. This additive structure means the Poisson term integrates into the existing eigendecomposition-based optimizer without restructuring the per-edge likelihood evaluation.

### Assumptions and their consequences

**Equal-weight events.** A 100-base deletion counts the same as a 1-base deletion. This means the model cannot distinguish between a branch with one long indel and a branch with one short indel. For typical viral phylodynamics (TreeTime's primary use case), indels are rare and short, so this has negligible impact. For bacterial genomics or gene family evolution with frequent long indels, the model underweights long events. See [proposal for length-weighted extensions](../../port-proposals/optimize-indel-model-alternatives.md).

**Symmetric ins/del.** Insertions and deletions share a single rate. Empirical data consistently shows deletion rate exceeds insertion rate (see Chapter 5, SIM/RIM section). The symmetric assumption is acceptable for the branch-length-prevents-zero use case because the Poisson optimum $t^* = k/\mu$ depends on total event count, not on the ins/del ratio.

**Constant rate across tree.** No branch-specific variation. The rate is a global average, re-estimated each round.

### Implementation

v1: [`packages/treetime/src/commands/optimize/optimize_indel.rs`](../../../packages/treetime/src/commands/optimize/optimize_indel.rs) (Poisson log-likelihood and rate estimation), [`packages/treetime/src/commands/optimize/optimize_unified.rs`](../../../packages/treetime/src/commands/optimize/optimize_unified.rs) (Newton integration).

### Related known issues

- [Grid search zero-comparison ignores indel likelihood](../../port-known-issues/M-optimize-grid-zero-ignores-indels.md)
- [Timetree branch length distribution ignores indels](../../port-known-issues/N-timetree-branch-length-distribution-ignores-indels.md)

### Cross-references

- [Intentional change](../../port-intentional-changes/optimize-indel-contribution-to-likelihood.md)
- [Alternatives proposal](../../port-proposals/optimize-indel-model-alternatives.md)

---

## Affine gap penalty

Assigns a cost per indel event (gap opening penalty $d$) plus a per-position extension cost ($e$), where $d > e$. The $O(MN)$ DP algorithm was introduced by <a id="cite-2"></a>[Gotoh 1982](<https://doi.org/10.1016/0022-2836(82)90398-9>) [[2](#ref-2)] for pairwise sequence alignment. Affine gap penalties are a scoring heuristic for alignment algorithms, not an evolutionary model: they define a cost function for placing gaps in an alignment, not a probability distribution over evolutionary events.

### Mathematical formulation

For a gap of length $g$: $\text{cost}(g) = d + (g - 1) \cdot e$

This cost function is equivalent to a log-probability under a <a id="gloss-use-2"></a>geometric distribution <sup>[2](#gloss-2)</sup> on gap lengths: $P(g) \propto \exp(-d) \cdot \exp(-e)^{g-1}$, with the ratio $e/d$ controlling the expected length. The alignment DP uses three matrices (match, insert-in-X, insert-in-Y) to track gap state transitions.

### Why affine gaps are not an evolutionary model

<a id="cite-3"></a>[Rivas 2005](https://doi.org/10.1186/1471-2105-6-63) [[3](#ref-3)] proved a fundamental limitation: single insertion events with geometric instantaneous length distributions do NOT produce geometric insert lengths at finite evolutionary times. The finite-time gap length distribution is a convolution of multiple overlapping indel events, which makes it heavier-tailed than geometric. No consistent continuous-time Markov process produces exactly affine gap costs. The affine model is an approximation that becomes worse at larger evolutionary distances where multiple overlapping indels are common.

### Empirical indel length distributions

The geometric distribution implied by affine gaps does not match observed data:

- <a id="cite-4"></a>[Qian and Goldstein 2001](https://doi.org/10.1002/prot.1129) [[4](#ref-4)]: structural alignment gap lengths are multi-exponential (four components), not geometric
- <a id="cite-5"></a>[Cartwright 2009](https://doi.org/10.1093/molbev/msn275) [[5](#ref-5)]: <a id="gloss-use-1"></a>Zipf <sup>[1](#gloss-1)</sup> (power-law) model fits far better than geometric across multiple organisms
- <a id="cite-6"></a>[Wygoda et al. 2024](https://doi.org/10.1093/bioinformatics/btae043) [[6](#ref-6)]: <a id="gloss-use-7"></a>ABC <sup>[7](#gloss-7)</sup> model-selection across multiple datasets confirms Zipf fits better than geometric in most cases. No efficient alignment programs exist for Zipf-distributed indel lengths.

### Role in phylogenetic software

Affine gap penalties are used in alignment tools (MAFFT, MUSCLE, ClustalW) to produce the fixed alignment that phylogenetic inference then operates on. They are not used in the tree likelihood computation. The gap penalties influence branch length estimates indirectly: different gap penalties produce different alignments, which produce different substitution patterns, which produce different branch lengths. But the gap penalties themselves do not enter the likelihood function.

PRANK is a special case: it uses a phylogeny-aware progressive alignment strategy that places gaps according to the tree topology, but its gap penalties are still heuristic, not derived from a birth-death model. PRANK defaults to inferring ancestral residues as absent when ambiguous, which reduces the tendency to over-align unrelated regions (a common problem with other progressive aligners).

---

## Glossary

1. <a id="gloss-1"></a> **Zipf distribution (power law).** A discrete probability distribution where $P(k) \propto k^{-a}$ for exponent $a > 1$. Empirical indel length data fits Zipf better than geometric in most datasets ([Cartwright 2009](https://doi.org/10.1093/molbev/msn275) [[5](#ref-5)]; [Wygoda et al. 2024](https://doi.org/10.1093/bioinformatics/btae043) [[6](#ref-6)]). [↩](#gloss-use-1)
2. <a id="gloss-2"></a> **Geometric distribution.** $P(k) = (1-p)p^{k-1}$. The distribution implied by affine gap penalties. [Rivas 2005](https://doi.org/10.1186/1471-2105-6-63) [[3](#ref-3)] proved that geometric instantaneous rates do not produce geometric finite-time distributions. [↩](#gloss-use-2)
3. <a id="gloss-7"></a> **ABC (Approximate Bayesian Computation).** Inference by simulation and summary statistic comparison, bypassing likelihood computation. [↩](#gloss-use-7)

## References

1. <a id="ref-1"></a> Warnow, Tandy. 2012. "Standard Maximum Likelihood Analyses of Alignments with Gaps Can Be Statistically Inconsistent." _PLOS Currents Tree of Life_. https://doi.org/10.1371/currents.rrn1308
2. <a id="ref-2"></a> Gotoh, Osamu. 1982. "An Improved Algorithm for Matching Biological Sequences." _Journal of Molecular Biology_ 162(3):705-708. https://doi.org/10.1016/0022-2836(82)90398-9
3. <a id="ref-3"></a> Rivas, Elena. 2005. "Evolutionary Models for Insertions and Deletions in a Probabilistic Modeling Framework." _BMC Bioinformatics_ 6:63. https://doi.org/10.1186/1471-2105-6-63
4. <a id="ref-4"></a> Qian, Bin, and Richard A. Goldstein. 2001. "Distribution of Indel Lengths." _Proteins_ 45(1):102-104. https://doi.org/10.1002/prot.1129
5. <a id="ref-5"></a> Cartwright, Reed A. 2009. "Problems and Solutions for Estimating Indel Rates and Length Distributions." _MBE_ 26(2):473-480. https://doi.org/10.1093/molbev/msn275
6. <a id="ref-6"></a> Wygoda, Elya, et al. 2024. "Statistical Framework to Determine Indel-Length Distribution." _Bioinformatics_ 40(2):btae043. https://doi.org/10.1093/bioinformatics/btae043
7. <a id="ref-27"></a> Truszkowski, Jakub, and Nick Goldman. 2015. "Maximum Likelihood Phylogenetic Inference is Consistent on Multiple Sequence Alignments, with or without Gaps." _Systematic Biology_ 65(2):328-333. https://doi.org/10.1093/sysbio/syv089

See [consolidated references](references.md) for the complete bibliography.
