# GTR Substitution Models

## Substitution Models

All models are continuous-time Markov chains on the nucleotide (or amino acid) alphabet. Each defines a rate matrix Q such that the transition probability matrix over branch length t is P(t) = exp(Qt). Models are normalized so that the expected rate of change at equilibrium equals 1: `beta = 1 / (-sum_i pi_i * Q_ii)`, and Q is scaled by beta so that branch length directly represents expected substitutions per site.

### Conventions

v1 uses the column convention throughout. The probability $p_i$ of state $i$ at one site evolves as

$$\frac{dp_i}{dt} = \sum_j Q_{ij} p_j$$

where $Q_{ij}$ is the rate from state $j$ to state $i$. Probability conservation requires every column to sum to zero, $\sum_i Q_{ij} = 0$, so the diagonal is $Q_{jj} = -\sum_{i \neq j} Q_{ij}$. The solution is $p(t) = e^{\mathbf{Q}t} p(0)$, so $P_{ij}(t) = \left(e^{\mathbf{Q}t}\right)_{ij}$ is the probability of state $i$ after time $t$ given state $j$ at the start. Each column of $P(t)$ sums to 1.

Time-reversible models satisfy detailed balance $Q_{ij}\pi_j = Q_{ji}\pi_i$: in equilibrium the flux from $j$ to $i$ equals the flux from $i$ to $j$. v1 enforces it by construction with a symmetric exchangeability matrix $W$:

$$Q_{ij} = \pi_i W_{ij} \quad (i \neq j)$$

The average rate is $-\sum_i \pi_i Q_{ii} = \sum_{ij} \pi_i W_{ij} \pi_j$. The constructor `GTR::new()` [`packages/treetime/src/gtr/gtr.rs#L31-L114`](../../packages/treetime/src/gtr/gtr.rs#L31-L114) symmetrizes $W$, zeroes its diagonal, normalizes $\pi$, divides $W$ by the average rate (`avg_transition()` [`packages/treetime/src/gtr/gtr.rs#L223-L225`](../../packages/treetime/src/gtr/gtr.rs#L223-L225)), and sets the scalar rate `mu` according to the rate scale below. The stored $W$ therefore has average rate 1, and `mu` carries the overall scale: $P(t) = e^{\mu \mathbf{Q} t}$.

The builder argument `mu` is a `Mu` ([`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)) that carries its unit, so the unit cannot differ from the value. With `Mu::FromMatrix(mu)`, `mu` multiplies the given $W$, so the model keeps the absolute rates of $W$; custom, inferred, and mugration models use it, as v0 `GTR.custom()` and `GTR.infer()` do through `assign_rates()`, and the inferred rate of `infer_gtr_from_counts()` is a `Mu::FromMatrix` at its producer. The named models (JC69, K80, F81, HKY85, T92, TN93, JTT92) use `Mu::Normalized(mu)`: `mu` is the average substitution rate of the finished model, 1 by default, so branch lengths are expected substitutions per site. v0 does the same in two steps: the model constructor folds the scale of $W$ into `mu`, and `GTR.standard()` then sets `mu` to the requested rate [`packages/legacy/treetime/treetime/gtr.py#L460`](../../packages/legacy/treetime/treetime/gtr.py#L460). The stored rate is one `f64` in the unit of the finished model.

A run names its sequence model with one `enum GtrModel` value ([`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)): `Infer`, or a named model with its parameter struct (`JC69Params`, `K80Params`, `F81Params`, `HKY85Params`, `T92Params`, `TN93Params`, `Jtt92Params`). The command layer builds it from `--model` and `--model-params` with `fn parse_model_params()` ([`packages/app-commands/src/commands/shared/model.rs`](../../packages/app-commands/src/commands/shared/model.rs)). It accepts the keys of the v0 model constructors, as v0 `create_gtr()` does [`packages/legacy/treetime/treetime/wrappers.py#L36-L79`](../../packages/legacy/treetime/treetime/wrappers.py#L36-L79): `mu` for every named model, `kappa` for K80, HKY85 and T92, `kappa1` and `kappa2` for TN93, `pi` (also `pis`, `Pi`, `Pis`) for F81, HKY85 and TN93, and `pi_GC` for T92. Where v0 skips a malformed item and passes an unknown key on to the constructor, which ignores it, v1 rejects a malformed or repeated item, an unknown key, a frequency list of other than 4 entries, and a value outside its domain (a non-positive rate, `pi_GC` outside (0, 1)). The partition builders call `GtrModel::build()` with the alphabet of the partition, which replaces the alphabet field of the parameter struct; `Infer` has no parameters and is fitted from the Fitch counts instead.

Consumers of $P(t)$:

- `GTR::propagate_profile()` computes `profile * P(t)`, the distribution of descendant states given a distribution of ancestral states (rows of `profile` are sites)
- `GTR::evolve()` computes `profile * P(t)^T`, the likelihood of each ancestral state given a likelihood vector of descendant states

Both are at [`packages/treetime/src/gtr/gtr.rs#L154-L192`](../../packages/treetime/src/gtr/gtr.rs#L154-L192). They choose the computation by the shape of the profile. A profile with fewer rows than states, such as the one-row profile of a discrete trait, is multiplied with the eigenvectors, scaled column by column with $e^{\mu \lambda_k t}$, and multiplied back (`fn in_eigenbasis()`, [`packages/treetime/src/gtr/gtr.rs#L259-L273`](../../packages/treetime/src/gtr/gtr.rs#L259-L273)), which costs O(k^2) per row; models with site rates always take this path. A one-row profile goes through vector-matrix products (BLAS `dgemv`), because OpenBLAS `dgemm` first copies ("packs") both operands into work buffers, and for one row this copy took about a third of the CPU time of mugration. A profile with more rows, such as a nucleotide alignment, is multiplied with $P(t)$ from `expQt()`, which costs O(k^3) per branch and then O(k^2) per row. The eigenbasis path clamps negative round-off of the result instead of the matrix, so the two paths agree to rounding. One helper, `rate_eigvals_at()` [`packages/treetime/src/gtr/gtr.rs#L198-L200`](../../packages/treetime/src/gtr/gtr.rs#L198-L200), forms the exponents $\mu r t \lambda_k$ (with a rate factor $r$, 1 without one) for `expQt()`, `expQt_with_rate()`, `evolve()`, `propagate_profile()`, and `rate_eigvals()`.

The model state is private and read through accessors. `GTR::new()` rejects a non-finite or negative `mu`, both as given and after the rate scale resolves it (`Mu::FromMatrix` multiplies it by the average rate, which can overflow), a `pi` with a zero, negative, or non-finite entry or a non-finite sum (the symmetric eigendecomposition divides by $\sqrt{\pi_i}$), a `W` with a negative or non-finite entry, a zero average rate, and site rates that are negative or non-finite. `set_mu()` is the only update; it validates `mu` like the constructor and keeps the eigendecomposition, because the stored eigenvalues belong to the normalized $W$ and `mu` scales them only in the exponent.

The model does not know the alignment, so the length of the site rates is checked where the model meets a partition: the backward and forward marginal passes of the dense, sparse, and discrete partitions call `GTR::check_site_count()` [`packages/treetime/src/gtr/gtr.rs#L121-L129`](../../packages/treetime/src/gtr/gtr.rs#L121-L129) and fail when the number of site rates differs from the number of sites (one site for a discrete partition). Without the check, a dense profile would broadcast a single site rate over all sites, and a sparse pass would index past the end of a short vector.

v0 uses the same convention. Phylogenetics textbooks usually use the transpose (row convention, $Q_{ij}$ is the rate from $i$ to $j$, rows sum to zero). Formulas taken from such sources must be transposed before they are compared with v1 code.

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

Simplest model: equal equilibrium frequencies and equal substitution rates. Zero free parameters after normalization. For an alphabet of $n$ states, $\pi_i = 1/n$ and the normalized rates are $Q_{ij} = 1/(n-1)$ for $i \neq j$ and $Q_{ii} = -1$. The eigenvalues are $0$ (once, with eigenvector $\pi$) and $-n/(n-1)$ ($n-1$ times, degenerate). The transition probability has the closed form

$$P_{ij}(t) = \frac{1 - e^{-nt/(n-1)}}{n} + e^{-nt/(n-1)}\,\delta_{ij}$$

For nucleotides ($n = 4$): `P_ii(t) = 1/4 + 3/4 * exp(-4t/3)`, `P_ij(t) = 1/4 - 1/4 * exp(-4t/3)`, eigenvalues {0, -4/3, -4/3, -4/3}.

The closed form would allow propagating a profile with $n$ numbers per site instead of a matrix-vector product. v1 does not use it: `jc69()` (`#jc69`) at [`packages/treetime/src/gtr/get_gtr.rs#L120-L132`](../../packages/treetime/src/gtr/get_gtr.rs#L120-L132) builds a generic `GTR` over the canonical states of the alphabet, and propagation goes through the numerical eigendecomposition like every other model. JC69 has no model-specific code: the single decay rate certificate below applies to it as to every model with equal exchangeabilities.

### K80 (<a id="cite-2"></a>[Kimura 1980](https://doi.org/10.1007/BF01731581) [[2](#ref-2)])

Distinguishes transitions (purine-purine A<->G, pyrimidine-pyrimidine C<->T) from transversions (purine-pyrimidine changes). One free parameter: the transition/transversion ratio kappa. Equal base frequencies. Closed-form P(t) with two exponential terms.

`k80()` (`#k80`) at [`packages/treetime/src/gtr/get_gtr.rs#L143-L155`](../../packages/treetime/src/gtr/get_gtr.rs#L143-L155).

### F81 (<a id="cite-3"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[3](#ref-3)])

Unequal equilibrium frequencies, equal exchangeabilities. Generalizes JC69 by allowing non-uniform base composition. Three free parameters (three independent frequencies; fourth constrained to sum to 1). In the column convention, `P_ij(t) = pi_i * (1 - exp(-beta*t))` for i != j, with `beta = 1 / (1 - sum_i pi_i^2)`.

`f81()` (`#f81`) at [`packages/treetime/src/gtr/get_gtr.rs#L169-L181`](../../packages/treetime/src/gtr/get_gtr.rs#L169-L181). Accepts optional `pi` parameter for non-uniform frequencies.

### HKY85 (<a id="cite-4"></a>[Hasegawa, Kishino, and Yano 1985](https://doi.org/10.1007/BF02101694) [[4](#ref-4)])

Combines K80's transition/transversion distinction with F81's unequal base frequencies. Four free parameters (kappa + three independent frequencies). Before normalization, `Q_ij = kappa * pi_i` for transitions and `pi_i` for transversions (column convention). Closed-form P(t) with three distinct exponential terms; eigenvalues involve kappa and purine/pyrimidine frequency sums (pi_R = pi_A + pi_G, pi_Y = pi_C + pi_T).

`hky85()` (`#hky85`) at [`packages/treetime/src/gtr/get_gtr.rs#L195-L214`](../../packages/treetime/src/gtr/get_gtr.rs#L195-L214). Accepts optional `pi` parameter.

### T92 (<a id="cite-9"></a>[Tamura 1992](https://doi.org/10.1093/oxfordjournals.molbev.a040059) [[9](#ref-9)])

GC-content parameterization: a simplification of HKY85 enforcing Chargaff's second parity rule (pi_A = pi_T, pi_C = pi_G). Parameterized by a single GC-content value theta = pi_G + pi_C, reducing three frequency parameters to one.

`t92()` (`#t92`) at [`packages/treetime/src/gtr/get_gtr.rs#L231-L254`](../../packages/treetime/src/gtr/get_gtr.rs#L231-L254).

### TN93 (<a id="cite-5"></a>[Tamura and Nei 1993](https://doi.org/10.1093/oxfordjournals.molbev.a040023) [[5](#ref-5)])

Distinguishes the two transition types: purine transitions (A<->G, rate kappa_1) and pyrimidine transitions (C<->T, rate kappa_2). Five free parameters. An analytical eigendecomposition exists.

`tn93()` (`#tn93`) at [`packages/treetime/src/gtr/get_gtr.rs#L346-L376`](../../packages/treetime/src/gtr/get_gtr.rs#L346-L376).

### JTT92 (<a id="cite-6"></a>[Jones, Taylor, and Thornton 1992](https://doi.org/10.1093/bioinformatics/8.3.275) [[6](#ref-6)])

Empirical 20x20 amino acid substitution matrix derived from a large database of protein sequence alignments. The exchangeability parameters and equilibrium frequencies are fixed to empirically observed values rather than estimated from data.

`jtt92()` (`#jtt92`) at [`packages/treetime/src/gtr/get_gtr.rs#L271-L325`](../../packages/treetime/src/gtr/get_gtr.rs#L271-L325).

---

## Matrix Exponentiation

Computing P(t) = exp(Qt) for the general time-reversible model (<a id="cite-7"></a>[Felsenstein 2003](https://openlibrary.org/books/OL8127124M/Inferring_Phylogenies) [[7](#ref-7)]; <a id="cite-8"></a>[Moler and Van Loan 2003](https://doi.org/10.1137/S00361445024180) [[8](#ref-8)]) requires eigendecomposition of the rate matrix. Analytical closed-form expressions exist for the simpler models (JC69 through TN93), but v1 builds every model, presets included, through the same numerical eigendecomposition.

### Symmetrization trick

Eigenvectors of a general matrix are numerically difficult. Detailed balance reduces the problem to a symmetric matrix. With the diagonal matrix $\mathbf{D} = \mathrm{diag}(\sqrt{\pi})$ and the column convention $Q_{ij} = \pi_i W_{ij}$:

$$\tilde{\mathbf{Q}} = \mathbf{D}^{-1} \mathbf{Q} \mathbf{D}, \qquad \tilde{Q}_{ij} = \sqrt{\pi_i}\, W_{ij} \sqrt{\pi_j} \quad (i \neq j), \qquad \tilde{Q}_{ii} = Q_{ii}$$

$\tilde{\mathbf{Q}}$ is real and symmetric, so `eigh` gives real eigenvalues $\lambda_k$ and orthonormal eigenvectors $w^k$. Because $\mathbf{Q}\mathbf{D}w^k = \mathbf{D}\tilde{\mathbf{Q}}w^k = \lambda_k \mathbf{D}w^k$, the right eigenvectors of $\mathbf{Q}$ are $\mathbf{D}w^k$ and the left eigenvectors are $\mathbf{D}^{-1}w^k$. The two sets are biorthonormal, and

$$P(t) = \mathbf{D}\, \mathbf{V} \,\mathrm{diag}\!\left(e^{\mu \lambda_k t}\right) \mathbf{V}^T \mathbf{D}^{-1}$$

where the columns of $\mathbf{V}$ are the $w^k$. In the row convention used by most textbooks the same matrix is written $\Pi^{1/2} Q \Pi^{-1/2}$.

`eig_single_site()` stores the right eigenvectors in `v`, scaled to unit L1 norm, and the left eigenvectors in `v_inv`, scaled by the inverse factor so that `v * v_inv = I`. `expQt()` computes `v * diag(exp(mu * lambda * t)) * v_inv` and clamps negative round-off to 0.

The eigendecomposition is computed once per model. Per branch, $P(t)$ needs only the diagonal exponentials and the stored eigenvector matrices, O(k^3) for the product instead of a matrix exponential; `propagate_profile()` and `evolve()` avoid the product for profiles with fewer rows than states (see above).

v1: `eig_single_site()` (`#eig_single_site`) at [`packages/treetime/src/gtr/gtr.rs#L273-L294`](../../packages/treetime/src/gtr/gtr.rs#L273-L294), `expQt()` (`#expQt`) at [`packages/treetime/src/gtr/gtr.rs#L194-L196`](../../packages/treetime/src/gtr/gtr.rs#L194-L196).

One eigenvalue is always 0. Its right eigenvector is the stationary distribution $\pi$, and its left eigenvector is the all-ones vector, because every column of $\mathbf{Q}$ sums to zero. Biorthogonality then implies that every other right eigenvector sums to zero. For an irreducible model the remaining k-1 eigenvalues are negative, guaranteeing convergence to equilibrium frequencies as t -> infinity.

### Single decay rate certificate

When all off-diagonal exchangeabilities are positive and equal, $W_{ij} = w$ for $i \neq j$ (the F81 family: binary models, JC69, F81, K80 and HKY85 with kappa 1, TN93 with its default `kappa1 = kappa2 = 1`, and custom, inferred, or mugration models with equal exchangeabilities), the rate matrix is $\mathbf{Q} = w(\pi \mathbf{1}^T - \mathbf{I})$ with eigenvalues $0$ and $-w$ ($k-1$ times), for any $\pi$ and any number of states. With the normalized $W$, $w = 1/(1 - \sum_i \pi_i^2)$. The branch log likelihood of an edge is then

$$\ell(t) = F(x) = \sum_s m_s \log(A_s + B_s x), \qquad x = e^{-a t}, \quad a = \mu w$$

$F$ is concave in $x$, and $x$ falls from 1 as $t$ grows, so $\ell'(0) < 0$ implies $\ell(t) \le \ell(0)$ for all $t \ge 0$: the optimal branch length is 0.

`GTR::new()` decides the certificate from the structure of $W$ after symmetrization and before the eigendecomposition, with exact equality and no tolerance [`packages/treetime/src/gtr/gtr.rs#L249-L254`](../../packages/treetime/src/gtr/gtr.rs#L249-L254). `certified_decay_rate()` returns $a$ for certified models and `None` otherwise [`packages/treetime/src/gtr/gtr.rs#L131-L133`](../../packages/treetime/src/gtr/gtr.rs#L131-L133). The optimizer takes the zero branch length shortcut only when all contributions of an edge share one certified decay rate [`packages/treetime/src/optimize/zero_boundary.rs#L120-L127`](../../packages/treetime/src/optimize/zero_boundary.rs#L120-L127). Contributions with different decay rates are not covered by the argument: $L(t) = (1 + 0.96 e^{-10t})(1 - 0.8 e^{-t})$, a product of two binary models, has $\ell'(0) \approx -0.898$ but $L(0) = 0.392 < L(0.5) \approx 0.518$.

The structural test certifies every model whose exchangeabilities are equal, whatever its name. TN93 with its default `kappa1 = kappa2 = 1` is therefore certified, while K80, HKY85, and T92 with their default kappa 0.1 are not; `--model-params` with `kappa=1` (or `kappa1` and `kappa2` other than 1) changes the certification. Certifying TN93 lets the optimizer take the zero branch length shortcut for TN93 edges. The effect on outputs is expected to be nil, but it was not measured.

---

## GTR Inference

### Iterative Coordinate Descent

Infers GTR model parameters (exchangeability matrix W, equilibrium frequencies pi, rate mu) from observed substitution patterns on the tree. The algorithm iterates: update W from transition counts, normalize, update pi from state occupancy, update mu from total rate. This coordinate descent converges to a local maximum of the likelihood.

The `MutationCounts` (`#MutationCounts`) struct holds the sufficient statistics: `nij` (directed substitution counts, `nij[i][j]` counts changes from state j to state i), `Ti` (time-in-state vector), and `root_state` (root state counts). Both sparse and dense inference paths produce `MutationCounts`, then share the same `infer_gtr_from_counts()` solver.

The updates are the approximate maximum-likelihood updates of <a id="cite-10"></a>[Puller et al. 2020](https://doi.org/10.1093/ve/veaa066) [[10](#ref-10)], restricted to one $\pi$ and one $\mu$ for all sites. The paper derives them for site-specific $\pi^a$ and $\mu^a$ with a shared $W$; see [kb/reports/auto-partitioning.md](../reports/auto-partitioning.md) for the site-specific form.

v1: `infer_gtr_from_counts()` (`#infer_gtr_from_counts`) at [`packages/treetime/src/gtr/infer_gtr/common.rs#L14-L87`](../../packages/treetime/src/gtr/infer_gtr/common.rs#L14-L87).
v0: [`packages/legacy/treetime/treetime/gtr.py#L491-L599`](../../packages/legacy/treetime/treetime/gtr.py#L491-L599).

### Fitch GTR Inference

Counts mutations from Fitch reconstruction: integer substitution counts for `nij`, branch-length-weighted composition for `Ti`, root composition from consensus sequence. Fast because Fitch reconstruction gives hard assignments (no probabilistic profiles to integrate over). Used by both dense and sparse initial GTR inference via `PartitionFitch::infer_gtr`.

`infer_gtr_fitch()` (`#infer_gtr_fitch`) and `get_mutation_counts_fitch()` at [`packages/treetime/src/partition/fitch/gtr_inference.rs`](../../packages/treetime/src/partition/fitch/gtr_inference.rs).

### Dense GTR Inference

Counts mutations from fractional expected counts derived from branch joint distributions. Requires two marginal reconstruction passes to populate profiles before GTR inference can run (the profiles are the input).

Key functions:

- `get_branch_mutation_matrix()` (`#get_branch_mutation_matrix`) at [`packages/treetime/src/gtr/infer_gtr/common.rs#L132-L161`](../../packages/treetime/src/gtr/infer_gtr/common.rs#L132-L161): computes posterior `P(child=i, parent=j | site)` from edge messages and transition matrix.
- `accumulate_mutation_counts()` (`#accumulate_mutation_counts`) at [`packages/treetime/src/gtr/infer_gtr/common.rs#L163-L194`](../../packages/treetime/src/gtr/infer_gtr/common.rs#L163-L194): sums `nij` and `Ti` from branch joint distributions.
- `count_transitions_dense()` (`#count_transitions_dense`) at [`packages/treetime/src/partition/marginal/shared/data.rs#L16-L57`](../../packages/treetime/src/partition/marginal/shared/data.rs#L16-L57): iterates edges to build `MutationCounts` with `SUPERTINY_NUMBER` floor on expQt and the effective branch length of each edge.

### GTR Refinement

Refits the GTR to expected transition counts of the marginal reconstruction, starting from an initial model. Only mugration refines, as v0 does; the sequence commands fit the GTR once, from the Fitch counts. `fn refine_gtr()` at [`packages/treetime/src/gtr/refinement.rs#L17-L83`](../../packages/treetime/src/gtr/refinement.rs#L17-L83) runs a marginal update under the start model, then the fits and rate steps, for any partition that implements `trait TransitionCounts` ([`packages/treetime/src/partition/marginal/shared/update.rs#L75`](../../packages/treetime/src/partition/marginal/shared/update.rs#L75)). The indexed partitions (dense sequences and discrete traits) implement it through their shared passes; the sparse partition does not count transitions. `struct GtrRefinement { iterations, inference, sampling_bias_correction }` ([`packages/treetime/src/gtr/refinement.rs#L86-L90`](../../packages/treetime/src/gtr/refinement.rs#L86-L90)) holds the iteration count, the `InferGtrOptions` of each fit, and the optional sampling-bias scale. Each fit runs `count_transitions()` and `infer_gtr_from_counts()`; a final marginal update under the refined model ends the run.

The order is the v0 order of `reconstruct_discrete_traits()`, with a rate step after each fit and the sampling-bias scale before the final update; the step list is in [kb/algo/mugration.md](mugration.md#gtr-model-construction). The rate step `fn optimize_gtr_rate()` ([`packages/treetime/src/gtr/refinement.rs#L126-L234`](../../packages/treetime/src/gtr/refinement.rs#L126-L234)) runs backward passes only and minimizes the negative log likelihood in $\sqrt{\mu}$ with `BrentBracketed` ([`packages/treetime/src/gtr/brent_bracketed.rs`](../../packages/treetime/src/gtr/brent_bracketed.rs)) on the range from $0.01\sqrt{\mu}$ to $100\sqrt{\mu}$. Without a bracket, `RateRangeEdge::KeepStart` keeps the start rate together with the messages of the evaluation at the start rate ([kb/decisions/mugration-rate-step-without-bracket-keeps-start-messages.md](../decisions/mugration-rate-step-without-bracket-keeps-start-messages.md)).

`RateRangeEdge::LowerCostEnd` ([`packages/treetime/src/gtr/refinement.rs#L237-L240`](../../packages/treetime/src/gtr/refinement.rs#L237-L240)) belongs to the relative-rate step of runs with several alignment files, `with_estimated_rate()` of the dense and sparse reconstructions ([`packages/treetime/src/partition/marginal/reconstruction.rs#L344-L359`](../../packages/treetime/src/partition/marginal/reconstruction.rs#L344-L359)): without a bracket it takes the end of the range with the lower cost ([kb/decisions/multi-alignment-partitions.md](../decisions/multi-alignment-partitions.md)).

---

## GTR Output

The model document (`--output-model`, `.model.json` by default) is `struct ModelOutput { alignments, models }` at [`packages/treetime/src/gtr/get_gtr.rs#L37-L46`](../../packages/treetime/src/gtr/get_gtr.rs#L37-L46). `alignments` lists the alignment files of the run (`struct AlignmentModel`: alignment number, path, number of sites, and the indel rate of the file); mugration leaves it out. `models` holds one entry per partition with a model, tagged `kind`: a `sequence` entry (`struct SequenceModel`: `partition`, `name`, `alignment` (the alignment number), `sites`, `clock_rate` in timetree, `states`) or a `trait` entry (`struct TraitModel`: `attribute`, `states`), each with the parameters of `struct GtrOutput` (model type, model name, mu, pi, W, number of states) flattened into it. `states` lists the characters of the states in the order of `pi` and `W`: the canonical characters of the partition's alphabet (`fn model_states()`, [`packages/treetime/src/gtr/get_gtr.rs#L25`](../../packages/treetime/src/gtr/get_gtr.rs#L25)) or the trait states. The commands build the document with `fn partition_models()` and write it with `fn write_gtr_output()` ([`packages/app-commands/src/commands/shared/gtr_output.rs`](../../packages/app-commands/src/commands/shared/gtr_output.rs)). Parameters are logged at info level via `log_gtr()` (`#log_gtr`).

---

## Unimplemented

See [unimplemented](unimplemented.md) for full details:

- Site-specific GTR ([kb/issues/M-gtr-site-specific-model-not-implemented.md](../issues/M-gtr-site-specific-model-not-implemented.md))
- Random GTR generation
- File-based GTR loading

Not tracked in [unimplemented](unimplemented.md), because v0 has neither:

- Closed-form propagation for JC69 (see the JC69 section above), which would replace the matrix-vector product per site with $n$ numbers
- Substitution models that change along the tree. v1 has one model per partition, and every edge of the partition uses it

---

## References

- <a id="ref-1"></a>Jukes, Thomas H., and Charles R. Cantor. 1969. "Evolution of Protein Molecules." In _Mammalian Protein Metabolism,_ vol. 3, edited by H. N. Munro, 21-132. Academic Press. https://doi.org/10.1016/B978-1-4832-3211-9.50009-7 [↩](#cite-1)
- <a id="ref-2"></a>Kimura, Motoo. 1980. "A Simple Method for Estimating Evolutionary Rates of Base Substitutions Through Comparative Studies of Nucleotide Sequences." _Journal of Molecular Evolution_ 16(2):111-120. https://doi.org/10.1007/BF01731581 [↩](#cite-2)
- <a id="ref-3"></a>Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17(6):368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-3)
- <a id="ref-4"></a>Hasegawa, Masami, Hirohisa Kishino, and Taka-aki Yano. 1985. "Dating of the Human-Ape Splitting by a Molecular Clock of Mitochondrial DNA." _Journal of Molecular Evolution_ 22(2):160-174. https://doi.org/10.1007/BF02101694 [↩](#cite-4)
- <a id="ref-5"></a>Tamura, Koichiro, and Masatoshi Nei. 1993. "Estimation of the Number of Nucleotide Substitutions in the Control Region of Mitochondrial DNA in Humans and Chimpanzees." _Molecular Biology and Evolution_ 10(3):512-526. https://doi.org/10.1093/oxfordjournals.molbev.a040023 [↩](#cite-5)
- <a id="ref-6"></a>Jones, David T., William R. Taylor, and Janet M. Thornton. 1992. "The Rapid Generation of Mutation Data Matrices from Protein Sequences." _Computer Applications in the Biosciences_ 8(3):275-282. https://doi.org/10.1093/bioinformatics/8.3.275 [↩](#cite-6)
- <a id="ref-7"></a>Felsenstein, Joseph. 2003. _Inferring Phylogenies._ Sinauer Associates. ISBN 978-0-87893-177-4. https://openlibrary.org/books/OL8127124M/Inferring_Phylogenies [↩](#cite-7)
- <a id="ref-8"></a>Moler, Cleve, and Charles Van Loan. 2003. "Nineteen Dubious Ways to Compute the Exponential of a Matrix, Twenty-Five Years Later." _SIAM Review_ 45(1):3-49. https://doi.org/10.1137/S00361445024180 [↩](#cite-8)
- <a id="ref-9"></a>Tamura, Koichiro. 1992. "Estimation of the Number of Nucleotide Substitutions When There Are Strong Transition-Transversion and G+C-Content Biases." _Molecular Biology and Evolution_ 9(4):678-687. https://doi.org/10.1093/oxfordjournals.molbev.a040059 [↩](#cite-9)
- <a id="ref-10"></a>Puller, Vadim, Pavel Sagulenko, and Richard A. Neher. 2020. "Efficient Inference, Potential, and Limitations of Site-Specific Substitution Models." _Virus Evolution_ 6(2):veaa066. https://doi.org/10.1093/ve/veaa066 [↩](#cite-10)

---

## File Index

| File                                                                                                                       | Algorithms                                                            |
| -------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)                                               | GTR core, eigendecomposition, `expQt()` (`#expQt`)                    |
| [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)                                       | JC69, K80, F81, HKY85, T92, TN93, JTT92, GTR output JSON              |
| [`packages/treetime/src/gtr/infer_gtr/common.rs`](../../packages/treetime/src/gtr/infer_gtr/common.rs)                     | `MutationCounts`, `InferGtrOptions`, `infer_gtr_from_counts()`        |
| [`packages/treetime/src/partition/fitch/gtr_inference.rs`](../../packages/treetime/src/partition/fitch/gtr_inference.rs)   | Fitch GTR inference from parsimony mutation counts (dense and sparse) |
| [`packages/treetime/src/partition/marginal/shared/data.rs`](../../packages/treetime/src/partition/marginal/shared/data.rs) | Dense GTR inference from branch joint distributions                   |
| [`packages/treetime/src/gtr/refinement.rs`](../../packages/treetime/src/gtr/refinement.rs)                                 | GTR refinement `refine_gtr()`, rate step `optimize_gtr_rate()`        |
| [`packages/treetime/src/gtr/brent_bracketed.rs`](../../packages/treetime/src/gtr/brent_bracketed.rs)                       | Bracketed Brent minimizer of the rate step                            |
