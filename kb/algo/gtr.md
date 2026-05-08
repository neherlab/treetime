# GTR Substitution Models

[Back to index](README.md)

## Substitution Models

All models are continuous-time Markov chains on the nucleotide (or amino acid) alphabet. Each defines a rate matrix Q such that the transition probability matrix over branch length t is P(t) = exp(Qt). Models are normalized so that the expected rate of change at equilibrium equals 1: `beta = 1 / (-sum_i pi_i * Q_ii)`, and Q is scaled by beta so that branch length directly represents expected substitutions per site.

The models form a nested hierarchy where simpler models are special cases of more complex ones:

```
JC69 (0 free params)
  |-- K80 (1: add kappa)
  |     |-- HKY85 (4: add unequal pi)
  |           |-- TN93 (5: split kappa into kappa_1, kappa_2)
  |                 |-- GTR (8: all 6 exchangeabilities free)
  |-- F81 (3: add unequal pi)
        |-- HKY85 (4: add kappa)
```

All models listed below are implemented in [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs).

### JC69 (<a id="cite-1"></a>[Jukes and Cantor 1969](https://doi.org/10.1016/B978-1-4832-3211-9.50009-7) [[1](#ref-1)])

Simplest model: equal base frequencies (pi = 1/4 each), equal substitution rates. Zero free parameters after normalization. Admits a closed-form transition probability: `P_ii(t) = 1/4 + 3/4 * exp(-4t/3)`, `P_ij(t) = 1/4 - 1/4 * exp(-4t/3)`. Eigenvalues: {0, -4/3, -4/3, -4/3}.

`jc69()` (`#jc69`) at [`packages/treetime/src/gtr/get_gtr.rs#L183-L194`](../../packages/treetime/src/gtr/get_gtr.rs#L183-L194).

### K80 (<a id="cite-2"></a>[Kimura 1980](https://doi.org/10.1007/BF01731581) [[2](#ref-2)])

Distinguishes transitions (purine-purine A<->G, pyrimidine-pyrimidine C<->T) from transversions (purine-pyrimidine changes). One free parameter: the transition/transversion ratio kappa. Equal base frequencies. Closed-form P(t) with two exponential terms.

`k80()` (`#k80`) at [`packages/treetime/src/gtr/get_gtr.rs#L219-L225`](../../packages/treetime/src/gtr/get_gtr.rs#L219-L225).

### F81 (<a id="cite-3"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[3](#ref-3)])

Unequal equilibrium frequencies, equal exchangeabilities. Generalizes JC69 by allowing non-uniform base composition. Three free parameters (three independent frequencies; fourth constrained to sum to 1). `P_ij(t) = pi_j * (1 - exp(-beta*t))` for i != j.

`f81()` (`#f81`) at [`packages/treetime/src/gtr/get_gtr.rs#L247-L257`](../../packages/treetime/src/gtr/get_gtr.rs#L247-L257). Accepts optional `pi` parameter for non-uniform frequencies.

### HKY85 (<a id="cite-4"></a>[Hasegawa, Kishino, and Yano 1985](https://doi.org/10.1007/BF02101694) [[4](#ref-4)])

Combines K80's transition/transversion distinction with F81's unequal base frequencies. Four free parameters (kappa + three independent frequencies). `Q_ij = kappa * pi_j` for transitions, `pi_j` for transversions. Closed-form P(t) with three distinct exponential terms; eigenvalues involve kappa and purine/pyrimidine frequency sums (pi_R = pi_A + pi_G, pi_Y = pi_C + pi_T).

`hky85()` (`#hky85`) at [`packages/treetime/src/gtr/get_gtr.rs#L285-L298`](../../packages/treetime/src/gtr/get_gtr.rs#L285-L298). Accepts optional `pi` parameter.

### T92 (<a id="cite-9"></a>[Tamura 1992](https://doi.org/10.1093/oxfordjournals.molbev.a040059) [[9](#ref-9)])

GC-content parameterization: a simplification of HKY85 enforcing Chargaff's second parity rule (pi_A = pi_T, pi_C = pi_G). Parameterized by a single GC-content value theta = pi_G + pi_C, reducing three frequency parameters to one.

`t92()` (`#t92`) at [`packages/treetime/src/gtr/get_gtr.rs#L325-L342`](../../packages/treetime/src/gtr/get_gtr.rs#L325-L342).

### TN93 (<a id="cite-5"></a>[Tamura and Nei 1993](https://doi.org/10.1093/oxfordjournals.molbev.a040023) [[5](#ref-5)])

Distinguishes the two transition types: purine transitions (A<->G, rate kappa_1) and pyrimidine transitions (C<->T, rate kappa_2). Five free parameters. Closed-form P(t) with analytical eigendecomposition.

`tn93()` (`#tn93`) at [`packages/treetime/src/gtr/get_gtr.rs#L459-L487`](../../packages/treetime/src/gtr/get_gtr.rs#L459-L487).

### JTT92 (<a id="cite-6"></a>[Jones, Taylor, and Thornton 1992](https://doi.org/10.1093/bioinformatics/8.3.275) [[6](#ref-6)])

Empirical 20x20 amino acid substitution matrix derived from a large database of protein sequence alignments. The exchangeability parameters and equilibrium frequencies are fixed to empirically observed values rather than estimated from data.

`jtt92()` (`#jtt92`) at [`packages/treetime/src/gtr/get_gtr.rs#L360-L418`](../../packages/treetime/src/gtr/get_gtr.rs#L360-L418).

---

## Matrix Exponentiation

Computing P(t) = exp(Qt) for the general time-reversible model (<a id="cite-7"></a>[Felsenstein 2003](https://doi.org/10.1007/978-0-387-21337-7) [[7](#ref-7)]; <a id="cite-8"></a>[Moler and Van Loan 2003](https://doi.org/10.1137/S0036144502418150) [[8](#ref-8)]) requires eigendecomposition of the rate matrix. For simpler models (JC69 through TN93), analytical closed-form expressions exist. For GTR and inferred models, numerical eigendecomposition is required.

### Symmetrization trick

Time-reversible rate matrices satisfy the detailed balance condition `pi_i * Q_ij = pi_j * Q_ji`. This allows symmetrization: `S = Pi^{1/2} * Q * Pi^{-1/2}` produces a real symmetric matrix with guaranteed real eigenvalues and orthogonal eigenvectors. Decompose `S = V * D * V^T` (columns of V are eigenvectors), then transform back:

```
P(t) = Pi^{-1/2} * V * diag(exp(lambda_i * t)) * V^T * Pi^{1/2}
```

The eigendecomposition is computed once per model. Per-branch computation reduces to multiplying diagonal exponentials by pre-computed eigenvector matrices - O(k^2) per branch rather than a full matrix exponential.

v1: `eig_single_site()` (`#eig_single_site`) at [`packages/treetime/src/gtr/gtr.rs#L71-L102`](../../packages/treetime/src/gtr/gtr.rs#L71-L102), `expQt()` (`#expQt`) at [`packages/treetime/src/gtr/gtr.rs#L442-L452`](../../packages/treetime/src/gtr/gtr.rs#L442-L452).

One eigenvalue is always 0 (corresponding to the stationary distribution pi). The remaining k-1 eigenvalues are negative, guaranteeing convergence to equilibrium frequencies as t -> infinity.

---

## GTR Inference

### Iterative Coordinate Descent

Infers GTR model parameters (exchangeability matrix W, equilibrium frequencies pi, rate mu) from observed substitution patterns on the tree. The algorithm iterates: update W from transition counts, normalize, update pi from state occupancy, update mu from total rate. This coordinate descent converges to a local maximum of the likelihood.

The `MutationCounts` (`#MutationCounts`) struct holds the sufficient statistics: `nij` (transition count matrix, symmetric) and `Ti` (time-in-state vector). Both sparse and dense inference paths produce `MutationCounts`, then share the same `infer_gtr_impl()` solver.

v1: `infer_gtr_impl()` (`#infer_gtr_impl`) at [`packages/treetime/src/gtr/infer_gtr/common.rs#L102-L163`](../../packages/treetime/src/gtr/infer_gtr/common.rs#L102-L163).
v0: [`packages/legacy/treetime/treetime/gtr.py#L491-L599`](../../packages/legacy/treetime/treetime/gtr.py#L491-L599).

### Fitch GTR Inference

Counts mutations from Fitch reconstruction: integer substitution counts for `nij`, branch-length-weighted composition for `Ti`, root composition from consensus sequence. Fast because Fitch reconstruction gives hard assignments (no probabilistic profiles to integrate over). Used by both dense and sparse initial GTR inference via `PartitionFitch::infer_gtr`.

`infer_gtr_fitch()` (`#infer_gtr_fitch`) at [`packages/treetime/src/gtr/infer_gtr/fitch.rs`](../../packages/treetime/src/gtr/infer_gtr/fitch.rs).

### Dense GTR Inference

Counts mutations from fractional expected counts derived from branch joint distributions. Requires two marginal reconstruction passes to populate profiles before GTR inference can run (the profiles are the input).

Key functions:

- `get_branch_mutation_matrix()` (`#get_branch_mutation_matrix`) at [`packages/treetime/src/gtr/infer_gtr/dense.rs#L51-L82`](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L51-L82): computes posterior `P(child=i, parent=j | site)` from edge messages and transition matrix.
- `accumulate_mutation_counts()` (`#accumulate_mutation_counts`) at [`packages/treetime/src/gtr/infer_gtr/dense.rs#L92-L130`](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L92-L130): sums `nij` and `Ti` from branch joint distributions.
- `get_mutation_counts_dense()` (`#get_mutation_counts_dense`) at [`packages/treetime/src/gtr/infer_gtr/dense.rs#L138-L210`](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L138-L210): iterates edges to build `MutationCounts` with `SUPERTINY_NUMBER` floor on expQt and branch length clamping.

---

## GTR Output

`write_gtr_json()` (`#write_gtr_json`) at [`packages/treetime/src/gtr/get_gtr.rs#L81-L94`](../../packages/treetime/src/gtr/get_gtr.rs#L81-L94) writes GTR model parameters (model type, model name, mu, pi, W) to JSON. Accepts an optional `qualifier` parameter: `None` writes `gtr.json`, `Some("sparse")` writes `gtr_sparse.json`, `Some("dense")` writes `gtr_dense.json`. Commands with a single partition pass `None`; the `optimize` command passes partition-type qualifiers to avoid overwriting when both sparse and dense partitions coexist. Parameters are logged at info level via `log_gtr()` (`#log_gtr`).

---

## Site-Specific GTR

Extension of the standard GTR where equilibrium frequencies vary per alignment site, requiring per-site eigendecomposition. The shared exchangeability matrix W captures relative rates between states. Per-site pi values produce per-site rate matrices Q_a = f(W, pi_a), each with its own eigendecomposition.

`GTRSiteSpecific` (`#GTRSiteSpecific`) at [`packages/treetime/src/gtr/gtr_site_specific.rs`](../../packages/treetime/src/gtr/gtr_site_specific.rs).

Core operations:

- Per-site eigendecomposition via `eig_single_site()` called per position
- `expQt(t)` returns Array3 [n_states, n_states, seq_len] with optional linear interpolation
- `propagate_profile()` and `evolve()` for backward/forward message passing
- `infer_gtr_site_specific_impl()` at [`packages/treetime/src/gtr/infer_gtr/site_specific.rs`](../../packages/treetime/src/gtr/infer_gtr/site_specific.rs) for iterative W/pi/mu inference

Not yet integrated into the partition system. See [unimplemented](unimplemented.md) for remaining work.

---

## Unimplemented

See [unimplemented](unimplemented.md) for full details:

- Site-specific GTR partition integration (core math implemented, wiring to partition types pending)
- Random GTR generation
- File-based GTR loading

---

## References

- <a id="ref-1"></a>Jukes, Thomas H., and Charles R. Cantor. 1969. "Evolution of Protein Molecules." In _Mammalian Protein Metabolism,_ vol. 3, edited by H. N. Munro, 21-132. Academic Press. https://doi.org/10.1016/B978-1-4832-3211-9.50009-7 [↩](#cite-1)
- <a id="ref-2"></a>Kimura, Motoo. 1980. "A Simple Method for Estimating Evolutionary Rates of Base Substitutions Through Comparative Studies of Nucleotide Sequences." _Journal of Molecular Evolution_ 16(2):111-120. https://doi.org/10.1007/BF01731581 [↩](#cite-2)
- <a id="ref-3"></a>Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-3)
- <a id="ref-4"></a>Hasegawa, Masami, Hirohisa Kishino, and Taka-aki Yano. 1985. "Dating of the Human-Ape Splitting by a Molecular Clock of Mitochondrial DNA." _Journal of Molecular Evolution_ 22(2):160-174. https://doi.org/10.1007/BF02101694 [↩](#cite-4)
- <a id="ref-5"></a>Tamura, Koichiro, and Masatoshi Nei. 1993. "Estimation of the Number of Nucleotide Substitutions in the Control Region of Mitochondrial DNA in Humans and Chimpanzees." _Molecular Biology and Evolution_ 10(3):512-526. https://doi.org/10.1093/oxfordjournals.molbev.a040023 [↩](#cite-5)
- <a id="ref-6"></a>Jones, David T., William R. Taylor, and Janet M. Thornton. 1992. "The Rapid Generation of Mutation Data Matrices from Protein Sequences." _Computer Applications in the Biosciences_ 8(3):275-282. https://doi.org/10.1093/bioinformatics/8.3.275 [↩](#cite-6)
- <a id="ref-7"></a>Felsenstein, Joseph. 2003. _Inferring Phylogenies._ Sinauer Associates. ISBN 978-0-87893-177-4. [↩](#cite-7)
- <a id="ref-8"></a>Moler, Cleve, and Charles Van Loan. 2003. "Nineteen Dubious Ways to Compute the Matrix Exponential, Twenty-Five Years Later." _SIAM Review_ 45(1):3-49. https://doi.org/10.1137/S0036144502418150 [↩](#cite-8)
- <a id="ref-9"></a>Tamura, Koichiro. 1992. "Estimation of the Number of Nucleotide Substitutions When There Are Strong Transition-Transversion and G+C-Content Biases." _Molecular Biology and Evolution_ 9(4):678-687. https://doi.org/10.1093/oxfordjournals.molbev.a040059 [↩](#cite-9)

---

## File Index

| File                                                                                                                 | Algorithms                                                            |
| -------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)                                         | GTR core, eigendecomposition, `expQt()` (`#expQt`)                    |
| [`packages/treetime/src/gtr/gtr_site_specific.rs`](../../packages/treetime/src/gtr/gtr_site_specific.rs)             | Site-specific GTR, per-site eigendecomposition, interpolation         |
| [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)                                 | JC69, K80, F81, HKY85, T92, TN93, JTT92, GTR output JSON              |
| [`packages/treetime/src/gtr/infer_gtr/common.rs`](../../packages/treetime/src/gtr/infer_gtr/common.rs)               | `MutationCounts`, `InferGtrOptions`, `infer_gtr_impl()`               |
| [`packages/treetime/src/gtr/infer_gtr/fitch.rs`](../../packages/treetime/src/gtr/infer_gtr/fitch.rs)                 | Fitch GTR inference from parsimony mutation counts (dense and sparse) |
| [`packages/treetime/src/gtr/infer_gtr/dense.rs`](../../packages/treetime/src/gtr/infer_gtr/dense.rs)                 | Dense GTR inference from branch joint distributions                   |
| [`packages/treetime/src/gtr/infer_gtr/site_specific.rs`](../../packages/treetime/src/gtr/infer_gtr/site_specific.rs) | Site-specific GTR inference from per-site mutation counts             |
