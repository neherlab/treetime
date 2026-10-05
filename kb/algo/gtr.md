# GTR Substitution Models

## Substitution Models

All models are continuous-time Markov chains on the nucleotide (or amino acid) alphabet. Each defines a rate matrix Q such that the transition probability matrix over branch length t is P(t) = exp(Qt). Models are normalized so that the expected rate of change at equilibrium equals 1: `beta = 1 / (-sum_i pi_i * Q_ii)`, and Q is scaled by beta so that branch length directly represents expected substitutions per site.

### Conventions

v1 uses the column convention throughout. The probability $p_i$ of state $i$ at one site evolves as

$$\frac{dp_i}{dt} = \sum_j Q_{ij} p_j$$

where $Q_{ij}$ is the rate from state $j$ to state $i$. Probability conservation requires every column to sum to zero, $\sum_i Q_{ij} = 0$, so the diagonal is $Q_{jj} = -\sum_{i \neq j} Q_{ij}$. The solution is $p(t) = e^{\mathbf{Q}t} p(0)$, so $P_{ij}(t) = \left(e^{\mathbf{Q}t}\right)_{ij}$ is the probability of state $i$ after time $t$ given state $j$ at the start. Each column of $P(t)$ sums to 1.

Time-reversible models satisfy detailed balance $Q_{ij}\pi_j = Q_{ji}\pi_i$: in equilibrium the flux from $j$ to $i$ equals the flux from $i$ to $j$. v1 enforces it by construction with a symmetric exchangeability matrix $W$:

$$Q_{ij} = \pi_i W_{ij} \quad (i \neq j)$$

The average rate is $-\sum_i \pi_i Q_{ii} = \sum_{ij} \pi_i W_{ij} \pi_j$. The constructor `GTR::new()` [`packages/treetime/src/gtr/gtr.rs#L54-L72`](../../packages/treetime/src/gtr/gtr.rs#L54-L72) symmetrizes $W$, zeroes its diagonal, normalizes $\pi$, divides $W$ by the average rate (`avg_transition()` [`packages/treetime/src/gtr/gtr.rs#L154-L156`](../../packages/treetime/src/gtr/gtr.rs#L154-L156)), and multiplies the scalar rate `mu` by it. The stored $W$ therefore has average rate 1, and `mu` carries the overall scale: $P(t) = e^{\mu \mathbf{Q} t}$.

Consumers of $P(t)$:

- `GTR::propagate_profile()` computes `profile · P(t)`, the distribution of descendant states given a distribution of ancestral states (rows of `profile` are sites)
- `GTR::evolve()` computes `profile · P(t)^T`, the likelihood of each ancestral state given a likelihood vector of descendant states

Both are at [`packages/treetime/src/gtr/gtr.rs#L103-L135`](../../packages/treetime/src/gtr/gtr.rs#L103-L135).

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

The closed form would allow propagating a profile with $n$ numbers per site instead of a matrix-vector product. v1 does not use it: `jc69()` (`#jc69`) at [`packages/treetime/src/gtr/get_gtr.rs#L118-L126`](../../packages/treetime/src/gtr/get_gtr.rs#L118-L126) builds a generic `GTR` over the canonical states of the alphabet, and propagation goes through the numerical eigendecomposition like every other model. The only JC69-specific behavior is the flag `unimodal_branch_likelihood`, which enables the zero-length branch shortcut of the optimizer.

### K80 (<a id="cite-2"></a>[Kimura 1980](https://doi.org/10.1007/BF01731581) [[2](#ref-2)])

Distinguishes transitions (purine-purine A<->G, pyrimidine-pyrimidine C<->T) from transversions (purine-pyrimidine changes). One free parameter: the transition/transversion ratio kappa. Equal base frequencies. Closed-form P(t) with two exponential terms.

`k80()` (`#k80`) at [`packages/treetime/src/gtr/get_gtr.rs#L141-L147`](../../packages/treetime/src/gtr/get_gtr.rs#L141-L147).

### F81 (<a id="cite-3"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[3](#ref-3)])

Unequal equilibrium frequencies, equal exchangeabilities. Generalizes JC69 by allowing non-uniform base composition. Three free parameters (three independent frequencies; fourth constrained to sum to 1). In the column convention, `P_ij(t) = pi_i * (1 - exp(-beta*t))` for i != j, with `beta = 1 / (1 - sum_i pi_i^2)`.

`f81()` (`#f81`) at [`packages/treetime/src/gtr/get_gtr.rs#L165-L173`](../../packages/treetime/src/gtr/get_gtr.rs#L165-L173). Accepts optional `pi` parameter for non-uniform frequencies.

### HKY85 (<a id="cite-4"></a>[Hasegawa, Kishino, and Yano 1985](https://doi.org/10.1007/BF02101694) [[4](#ref-4)])

Combines K80's transition/transversion distinction with F81's unequal base frequencies. Four free parameters (kappa + three independent frequencies). Before normalization, `Q_ij = kappa * pi_i` for transitions and `pi_i` for transversions (column convention). Closed-form P(t) with three distinct exponential terms; eigenvalues involve kappa and purine/pyrimidine frequency sums (pi_R = pi_A + pi_G, pi_Y = pi_C + pi_T).

`hky85()` (`#hky85`) at [`packages/treetime/src/gtr/get_gtr.rs#L191-L204`](../../packages/treetime/src/gtr/get_gtr.rs#L191-L204). Accepts optional `pi` parameter.

### T92 (<a id="cite-9"></a>[Tamura 1992](https://doi.org/10.1093/oxfordjournals.molbev.a040059) [[9](#ref-9)])

GC-content parameterization: a simplification of HKY85 enforcing Chargaff's second parity rule (pi_A = pi_T, pi_C = pi_G). Parameterized by a single GC-content value theta = pi_G + pi_C, reducing three frequency parameters to one.

`t92()` (`#t92`) at [`packages/treetime/src/gtr/get_gtr.rs#L221-L238`](../../packages/treetime/src/gtr/get_gtr.rs#L221-L238).

### TN93 (<a id="cite-5"></a>[Tamura and Nei 1993](https://doi.org/10.1093/oxfordjournals.molbev.a040023) [[5](#ref-5)])

Distinguishes the two transition types: purine transitions (A<->G, rate kappa_1) and pyrimidine transitions (C<->T, rate kappa_2). Five free parameters. An analytical eigendecomposition exists.

`tn93()` (`#tn93`) at [`packages/treetime/src/gtr/get_gtr.rs#L328-L352`](../../packages/treetime/src/gtr/get_gtr.rs#L328-L352).

### JTT92 (<a id="cite-6"></a>[Jones, Taylor, and Thornton 1992](https://doi.org/10.1093/bioinformatics/8.3.275) [[6](#ref-6)])

Empirical 20x20 amino acid substitution matrix derived from a large database of protein sequence alignments. The exchangeability parameters and equilibrium frequencies are fixed to empirically observed values rather than estimated from data.

`jtt92()` (`#jtt92`) at [`packages/treetime/src/gtr/get_gtr.rs#L255-L303`](../../packages/treetime/src/gtr/get_gtr.rs#L255-L303).

---

## Matrix Exponentiation

Computing P(t) = exp(Qt) for the general time-reversible model (<a id="cite-7"></a>[Felsenstein 2003](https://doi.org/10.1007/978-0-387-21337-7) [[7](#ref-7)]; <a id="cite-8"></a>[Moler and Van Loan 2003](https://doi.org/10.1137/S0036144502418150) [[8](#ref-8)]) requires eigendecomposition of the rate matrix. Analytical closed-form expressions exist for the simpler models (JC69 through TN93), but v1 builds every model, presets included, through the same numerical eigendecomposition.

### Symmetrization trick

Eigenvectors of a general matrix are numerically difficult. Detailed balance reduces the problem to a symmetric matrix. With the diagonal matrix $\mathbf{D} = \mathrm{diag}(\sqrt{\pi})$ and the column convention $Q_{ij} = \pi_i W_{ij}$:

$$\tilde{\mathbf{Q}} = \mathbf{D}^{-1} \mathbf{Q} \mathbf{D}, \qquad \tilde{Q}_{ij} = \sqrt{\pi_i}\, W_{ij} \sqrt{\pi_j} \quad (i \neq j), \qquad \tilde{Q}_{ii} = Q_{ii}$$

$\tilde{\mathbf{Q}}$ is real and symmetric, so `eigh` gives real eigenvalues $\lambda_k$ and orthonormal eigenvectors $w^k$. Because $\mathbf{Q}\mathbf{D}w^k = \mathbf{D}\tilde{\mathbf{Q}}w^k = \lambda_k \mathbf{D}w^k$, the right eigenvectors of $\mathbf{Q}$ are $\mathbf{D}w^k$ and the left eigenvectors are $\mathbf{D}^{-1}w^k$. The two sets are biorthonormal, and

$$P(t) = \mathbf{D}\, \mathbf{V} \,\mathrm{diag}\!\left(e^{\mu \lambda_k t}\right) \mathbf{V}^T \mathbf{D}^{-1}$$

where the columns of $\mathbf{V}$ are the $w^k$. In the row convention used by most textbooks the same matrix is written $\Pi^{1/2} Q \Pi^{-1/2}$.

`eig_single_site()` stores the right eigenvectors in `v`, scaled to unit L1 norm, and the left eigenvectors in `v_inv`, scaled by the inverse factor so that `v · v_inv = I`. `expQt()` computes `v · diag(exp(mu * lambda * t)) · v_inv` and clamps negative round-off to 0.

The eigendecomposition is computed once per model. Per-branch computation reduces to multiplying diagonal exponentials by pre-computed eigenvector matrices - O(k^2) per branch rather than a full matrix exponential.

v1: `eig_single_site()` (`#eig_single_site`) at [`packages/treetime/src/gtr/gtr.rs#L158-L179`](../../packages/treetime/src/gtr/gtr.rs#L158-L179), `expQt()` (`#expQt`) at [`packages/treetime/src/gtr/gtr.rs#L137-L143`](../../packages/treetime/src/gtr/gtr.rs#L137-L143).

One eigenvalue is always 0. Its right eigenvector is the stationary distribution $\pi$, and its left eigenvector is the all-ones vector, because every column of $\mathbf{Q}$ sums to zero. Biorthogonality then implies that every other right eigenvector sums to zero. For an irreducible model the remaining k-1 eigenvalues are negative, guaranteeing convergence to equilibrium frequencies as t -> infinity.

---

## GTR Inference

### Iterative Coordinate Descent

Infers GTR model parameters (exchangeability matrix W, equilibrium frequencies pi, rate mu) from observed substitution patterns on the tree. The algorithm iterates: update W from transition counts, normalize, update pi from state occupancy, update mu from total rate. This coordinate descent converges to a local maximum of the likelihood.

The `MutationCounts` (`#MutationCounts`) struct holds the sufficient statistics: `nij` (directed substitution counts, `nij[i][j]` counts changes from state j to state i), `Ti` (time-in-state vector), and `root_state` (root state counts). Both sparse and dense inference paths produce `MutationCounts`, then share the same `infer_gtr_impl()` solver.

The updates are the approximate maximum-likelihood updates of <a id="cite-10"></a>[Puller et al. 2020](https://doi.org/10.1093/ve/veaa066) [[10](#ref-10)], restricted to one $\pi$ and one $\mu$ for all sites. The paper derives them for site-specific $\pi^a$ and $\mu^a$ with a shared $W$; see [kb/reports/auto-partitioning.md](../reports/auto-partitioning.md) for the site-specific form.

v1: `infer_gtr_impl()` (`#infer_gtr_impl`) at [`packages/treetime/src/gtr/infer_gtr/common.rs#L13-L82`](../../packages/treetime/src/gtr/infer_gtr/common.rs#L13-L82).
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

---

## GTR Output

`struct GtrOutput` at [`packages/treetime/src/gtr/get_gtr.rs#L27`](../../packages/treetime/src/gtr/get_gtr.rs#L27) holds the GTR model parameters (model type, model name, mu, pi, W). Each command writes it with `json_write_file()` to the path of the `gtr` output (`.gtr.json` by default). Parameters are logged at info level via `log_gtr()` (`#log_gtr`).

---

## Unimplemented

See [unimplemented](unimplemented.md) for full details:

- Site-specific GTR partition integration (core math implemented, wiring to partition types pending)
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
- <a id="ref-7"></a>Felsenstein, Joseph. 2003. _Inferring Phylogenies._ Sinauer Associates. ISBN 978-0-87893-177-4. [↩](#cite-7)
- <a id="ref-8"></a>Moler, Cleve, and Charles Van Loan. 2003. "Nineteen Dubious Ways to Compute the Matrix Exponential, Twenty-Five Years Later." _SIAM Review_ 45(1):3-49. https://doi.org/10.1137/S0036144502418150 [↩](#cite-8)
- <a id="ref-9"></a>Tamura, Koichiro. 1992. "Estimation of the Number of Nucleotide Substitutions When There Are Strong Transition-Transversion and G+C-Content Biases." _Molecular Biology and Evolution_ 9(4):678-687. https://doi.org/10.1093/oxfordjournals.molbev.a040059 [↩](#cite-9)
- <a id="ref-10"></a>Puller, Vadim, Pavel Sagulenko, and Richard A. Neher. 2020. "Efficient Inference, Potential, and Limitations of Site-Specific Substitution Models." _Virus Evolution_ 6(2):veaa066. https://doi.org/10.1093/ve/veaa066 [↩](#cite-10)

---

## File Index

| File                                                                                                                       | Algorithms                                                            |
| -------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| [`packages/treetime/src/gtr/gtr.rs`](../../packages/treetime/src/gtr/gtr.rs)                                               | GTR core, eigendecomposition, `expQt()` (`#expQt`)                    |
| [`packages/treetime/src/gtr/get_gtr.rs`](../../packages/treetime/src/gtr/get_gtr.rs)                                       | JC69, K80, F81, HKY85, T92, TN93, JTT92, GTR output JSON              |
| [`packages/treetime/src/gtr/infer_gtr/common.rs`](../../packages/treetime/src/gtr/infer_gtr/common.rs)                     | `MutationCounts`, `InferGtrOptions`, `infer_gtr_impl()`               |
| [`packages/treetime/src/partition/fitch/gtr_inference.rs`](../../packages/treetime/src/partition/fitch/gtr_inference.rs)   | Fitch GTR inference from parsimony mutation counts (dense and sparse) |
| [`packages/treetime/src/partition/marginal/shared/data.rs`](../../packages/treetime/src/partition/marginal/shared/data.rs) | Dense GTR inference from branch joint distributions                   |
