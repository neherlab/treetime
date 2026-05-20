# Chapter 2: Substitution models and sequence evolution

[Back to index](README.md) | Previous: [Chapter 1: Introduction](1-introduction.md) | Next: [Chapter 3: Tree likelihood](3-tree-likelihood.md)

## How DNA changes over time

DNA sequences evolve by **substitution**: one nucleotide replaces another at a given position. At position 42 in an alignment, an ancestor had nucleotide A; its descendant has nucleotide G. We write this as `A->G at position 42`.

Over evolutionary time, substitutions accumulate. The longer the time since two sequences shared a common ancestor, the more positions differ. This relationship is the basis of the molecular clock and of branch length estimation.

Not all substitutions are equally likely. Two empirical patterns dominate:

1. Transition/transversion bias. Transitions (purine-to-purine: A<->G, or pyrimidine-to-pyrimidine: C<->T) are 2x to 30x more frequent than transversions (purine-to-pyrimidine or vice versa). This reflects the biochemistry of DNA: transitions require only a single tautomeric shift, while transversions require a purine-pyrimidine swap.

2. Base composition bias. The four nucleotides are not equally frequent. AT-rich genomes (many bacteria, mitochondria) have more A and T; GC-rich genomes (some vertebrates) have more G and C. The equilibrium frequencies affect which substitutions are observed.

Substitution models formalize these patterns mathematically.

## The continuous-time Markov chain

A substitution model is a continuous-time Markov chain (CTMC) on four states {A, C, G, T}. The model is defined by a 4x4 **rate matrix** `Q`, where (row-stochastic convention: `Q_{ij}` = rate from i to j, rows sum to zero):

- `Q_{ij}` (off-diagonal): the instantaneous rate of substitution from state i to state j
- `Q_{ii}` (diagonal): set so each row sums to zero (`Q_{ii} = -sum_{j!=i} Q_{ij}`)

Note: the algorithm specification in `../_raw/sequence_evolution.md` uses the transpose (column-stochastic) convention where `Q_{ij}` is the rate from j to i and columns sum to zero. The two conventions are mathematically equivalent; this chapter uses the row-stochastic form common in phylogenetics textbooks.

The **transition probability matrix** at evolutionary distance t is the matrix exponential:

```
P(t) = exp(Q * t)
```

Entry `P(t)_{ij}` is the probability that state i at the ancestor becomes state j at the descendant after evolutionary distance t. At `t = 0`, `P(0) = I` (identity -- no change). As `t` grows, `P(t)` approaches the equilibrium distribution (all rows identical).

The Markov property means: the probability of the next state depends only on the current state, not on the history of how the current state was reached. This is the mathematical justification for treating each branch independently when computing tree likelihoods.

## The semigroup property and substitution composition

The transition matrices satisfy the Chapman-Kolmogorov equation:

```
P(t1 + t2) = P(t1) * P(t2)
```

This is **matrix multiplication**, not element-wise multiplication. The probability of going from state i to state j over time `t1 + t2` equals the sum over all intermediate states k:

```
P(t1+t2)_{ij} = sum_k P(t1)_{ik} * P(t2)_{kj}
```

This property has two critical consequences for tree refinement:

**Computing tree likelihoods.** The pruning algorithm ([Chapter 3](3-tree-likelihood.md)) multiplies transition matrices along branches from leaves to root.

**Collapsing edges.** When an internal node is removed (e.g., pruning a zero-length branch, [Chapter 6](6-zero-length-branches.md)), the substitutions on the two adjacent edges must be **composed**, not unioned. A concrete example:

```
Before collapse:             After collapse:

  grandparent (A at pos 5)     grandparent (A at pos 5)
       |                            |
   edge1: A->G at pos 5        single edge: A->T at pos 5
       |                            |
    parent (G at pos 5)         child (T at pos 5)
       |
   edge2: G->T at pos 5
       |
    child (T at pos 5)
```

Three cases arise when composing substitutions at the same position:

| Edge 1 | Edge 2 | Composition         | Set-union (wrong)             |
| ------ | ------ | ------------------- | ----------------------------- |
| `A->G` | `G->T` | `A->T` (chain)      | `{A->G, G->T}` (two entries)  |
| `A->G` | `G->A` | no change (cancel)  | `{A->G, G->A}` (two entries)  |
| `A->G` | `A->G` | `A->G` (convergent) | `A->G` (accidentally correct) |

Set-union produces correct results only for convergent mutations at the same position. For chains and cancellations, it produces incorrect mutation counts and incorrect ancestral state assignments.

v1 implements composition via `compose_substitutions()` in [`packages/treetime/src/seq/mutation.rs`](../../../packages/treetime/src/seq/mutation.rs), called from the shared `collapse_edge()` in [`packages/treetime/src/optimize/topology/collapse.rs`](../../../packages/treetime/src/optimize/topology/collapse.rs). Both the prune and optimize commands delegate to this shared function so sparse edge collapse stays composition-correct across the codebase.

## The model hierarchy

Substitution models form a hierarchy from simple to complex. Each model is obtained by relaxing constraints on the rate matrix Q. All are special cases of the GTR (General Time Reversible) model.

### JC69 (Jukes-Cantor 1969)

The simplest model. All substitution rates are equal, all base frequencies are 1/4. One free parameter: the overall rate mu.

```
Rate matrix Q (JC69):

       A      C      G      T
  A [ -3mu    mu     mu     mu  ]
  C [  mu   -3mu     mu     mu  ]
  G [  mu     mu   -3mu     mu  ]
  T [  mu     mu     mu   -3mu  ]
```

The transition probability has a closed form:

```
P(t)_{ii} = 1/4 + 3/4 * exp(-4*mu*t)     (same state)
P(t)_{ij} = 1/4 - 1/4 * exp(-4*mu*t)     (different state, i != j)
```

The per-edge log-likelihood is guaranteed unimodal under JC69: any optimization method converges to the global optimum (Dinh and Matsen 2017).

v1 code: JC69 is the default model. Constructed by `jc69()` in [`packages/treetime/src/gtr/get_gtr.rs`](../../../packages/treetime/src/gtr/get_gtr.rs).

### K2P (Kimura 1980)

Distinguishes transitions (rate alpha) from transversions (rate beta). Equal base frequencies. Two free parameters.

The transition/transversion ratio kappa = alpha/beta. For real data, kappa ranges from ~2 (mitochondrial DNA) to ~15 (some nuclear genes).

Under K2P, the per-edge likelihood can have **multiple local maxima**. Dinh and Matsen (2017) proved that the space of rescaled 1D likelihood functions under K2P is dense in all continuous non-negative functions on `[0, infinity)` -- the likelihood surface can take any shape. This means Newton-Raphson from a single starting point can converge to a local optimum.

### HKY85 (Hasegawa, Kishino, and Yano 1985)

Combines transition/transversion distinction with unequal base frequencies. Five free parameters: kappa + four frequencies (three free, summing to 1). The standard intermediate-complexity model, capturing both dominant empirical patterns.

### GTR (General Time Reversible)

The most general reversible model. The rate matrix factors as:

```
Q = S * diag(pi)
```

where `S` is a symmetric 6-parameter exchangeability matrix and `pi` is the stationary frequency vector. Eight free parameters total (6 exchangeabilities + 4 frequencies - 1 normalization - 1 frequency constraint).

All simpler models are special cases: JC69 constrains all exchangeabilities equal and all frequencies to 1/4; K2P allows two exchangeability classes; HKY allows two classes with unequal frequencies.

GTR was formulated independently by Lanave et al. (1984) and given its mathematical framework by Tavaré (1986).

v1 code: GTR inference in [`packages/treetime/src/gtr/get_gtr.rs`](../../../packages/treetime/src/gtr/get_gtr.rs). The `Gtr` struct stores the rate matrix, eigenvectors, eigenvalues, and stationary frequencies.

## The eigendecomposition

For any diagonalizable rate matrix Q, the eigendecomposition `Q = V * diag(lambda) * V^-1` enables efficient computation:

```
P(t) = V * diag(exp(lambda * t)) * V^-1
```

The per-site conditional likelihood for one edge becomes:

```
L_i(t) = sum_c k_{ic} * exp(lambda_c * t)
```

where the coefficients `k_{ic} = (msg_parent.dot(V))_c * (msg_child.dot(V_inv.T))_c` depend only on the ancestral state distributions at the endpoints, not on the branch length t. These coefficients are computed once per edge and reused across all Newton-Raphson iterations during branch length optimization ([Chapter 5](5-branch-length-optimization.md)).

The analytical first and second derivatives follow by multiplying by eigenvalues:

```
L'_i(t) = sum_c k_{ic} * lambda_c * exp(lambda_c * t)
L''_i(t) = sum_c k_{ic} * lambda_c^2 * exp(lambda_c * t)
```

This makes Newton-Raphson optimization efficient: function value and derivatives computed from the same precomputed coefficients, O(s) per site (s = number of states, 4 for nucleotides).

v1 code: eigendecomposition-based evaluation in [`packages/treetime/src/optimize/dense_eval.rs`](../../../packages/treetime/src/optimize/dense_eval.rs) (dense) and [`packages/treetime/src/optimize/sparse_eval.rs`](../../../packages/treetime/src/optimize/sparse_eval.rs) (sparse).

### Eigenvalues and model properties

For reversible models, all eigenvalues are non-positive. One eigenvalue is always zero (corresponding to the stationary distribution). The remaining eigenvalues determine how fast the Markov chain approaches equilibrium.

For JC69: eigenvalues are `[0, -4mu, -4mu, -4mu]` where `mu` is the per-pair substitution rate (the off-diagonal entry of Q), not the overall substitution rate. Three equal negative eigenvalues.

For GTR: eigenvalues are `[0, lambda_1, lambda_2, lambda_3]` with `lambda_i < 0`. Distinct for most parameter choices, making the likelihood surface more complex.

## Among-site rate variation

Not all alignment positions evolve at the same rate. Conserved positions (functional constraints) evolve slowly; variable positions (neutral sites) evolve faster. Yang (1994) introduced the discrete-gamma approximation: the continuous gamma distribution of site rates is approximated by K discrete categories (K=4) of equal probability. The site likelihood becomes:

```
L_site = (1/K) * sum_{k=1}^{K} L_site(r_k)
```

where `r_k` are the discrete rates. The shape parameter alpha controls rate variation: small alpha = high variation, alpha -> infinity = uniform rates.

This "+Gamma" extension can be combined with any substitution model (e.g., GTR+G). Branch length optimization must account for rate categories by summing over them.

v1 code: per-site rate variation is tracked as a known issue (`M-gtr-per-site-rate-variation`). Not yet implemented.

## References

- Jukes, T. H., and C. R. Cantor. 1969. "Evolution of Protein Molecules." In _Mammalian Protein Metabolism_, vol. 3, 21-132. https://doi.org/10.1016/B978-1-4832-3211-9.50009-7
- Kimura, M. 1980. "A Simple Method for Estimating Evolutionary Rates." _J. Mol. Evol._ 16(2):111-120. https://doi.org/10.1007/BF01731581
- Hasegawa, M., H. Kishino, and T. Yano. 1985. "Dating of the Human-Ape Splitting." _J. Mol. Evol._ 22(2):160-174. https://doi.org/10.1007/BF02101694
- Lanave, C., G. Preparata, C. Saccone, and G. Serio. 1984. "A New Method for Calculating Evolutionary Substitution Rates." _J. Mol. Evol._ 20:86-93. https://doi.org/10.1007/BF02101990
- Tavaré, S. 1986. "Some Probabilistic and Statistical Problems in the Analysis of DNA Sequences." _Lect. Math. Life Sci._ 17:57-86. No DOI available (pre-DOI publication).
- Yang, Z. 1994. "Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates." _J. Mol. Evol._ 39:306-314. https://doi.org/10.1007/BF00160154
- Dinh, V. C., and F. A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Ann. Appl. Prob._ 27(3):1646-1677. https://doi.org/10.1214/16-AAP1240
