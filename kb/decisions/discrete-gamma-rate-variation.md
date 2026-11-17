# Discrete gamma rate variation: K rate categories per partition with an estimated shape

`--gamma-categories K` (config key `gamma_categories`, K of 2 or more) gives every sequence partition of ancestral (marginal), optimize and timetree K rate categories of equal probability. The rates are the means of K equal-probability slices of a gamma distribution with mean 1 and shape α (Yang 1994), and each partition estimates α by maximum likelihood. The likelihood of a site is the mean of its likelihoods under the K rates, and the ancestral states of a site are weighted by the posterior probability of each rate at that site. Without the setting, all sites of a partition share one rate, as in v0.

**Type**: Scientific choice for v1 (v0 has no rate variation across sites).

**Status**: Approved on 2026-11-17.

## Rule

### Category rates

- **Mean-rate discretization**: category k of K covers the slice between the (k-1)/K and k/K quantiles a_k and b_k of the gamma distribution with shape α and rate α (mean 1). Its rate is the mean of the distribution over the slice, r_k = K * (P(α+1, α*b_k) - P(α+1, α*a_k)), with P the regularized lower incomplete gamma function. Every category has the weight 1/K (`struct GammaRates`, [packages/treetime/src/gtr/gamma.rs#L18-L58](../../packages/treetime/src/gtr/gamma.rs#L18-L58)). This is the default of RAxML-NG and IQ-TREE ([raxml-ng src/Model.cpp#L466](https://github.com/amkozlov/raxml-ng/blob/d396351ee1263b704adb135edd6a4a84522dbd5b/src/Model.cpp#L466), `.repos/iqtree3/model/rategamma.h`)
- **Mean 1**: the rates have mean 1, so branch lengths keep their unit of expected substitutions per site, and category k of partition p propagates with the scalar rate μ_p * r_k
- **Shape range**: α lies between 0.02 and 100, the bounds of coraxlib (`.repos/coraxlib/src/corax/optimize/opt_generic.h#L83-L84`); IQ-TREE uses 0.02 to 1000
- **Partition model**: the model of a partition is the pair of its GTR and its gamma rates (`struct PartitionModel`, [packages/treetime/src/gtr/partition_model.rs#L8](../../packages/treetime/src/gtr/partition_model.rs#L8)). The amino-acid reconstructions of `--translations` run without rate variation

### Setting

- One value for the whole run, applied to every partition; `--partition` entries have no gamma key
- It applies to ancestral with `--method-anc marginal`, and to optimize and timetree with `--alignment`. With `--method-anc parsimony` or `--branch-length-mode input` it is an argument error, because those runs build no marginal sequence partition
- Partitions share α only when they share `gamma-shape` (`--shared-params`, section "Shared model parameters" of [kb/decisions/multi-alignment-partitions.md](multi-alignment-partitions.md#shared-model-parameters))

### Reconstruction: one pass per category

- **Category passes**: a partition reconstruction with K >= 2 holds K passes, each the existing dense or sparse pass with the scalar rate μ_p * r_k, so the propagation code is the same as without gamma (`packages/treetime/src/partition/marginal/mixture/`). Sparse fixed positions are exact in every pass, because each pass uses one scalar rate. A run with K = 1 runs the plain pass
- **Per-site log-likelihoods**: each pass gives the root log-likelihood ℓ_{s,k} of every site s. Dense passes compute it after the backward pass from the stored child messages. Sparse passes keep a scale per state group and an offset per site during the backward pass; a fixed-state group counts once per site of the group
- **Mixture log-likelihood**: Σ_s log Σ_k (1/K) e^{ℓ_{s,k}}, with a log-sum-exp per site. It replaces the root log-likelihood wherever the reconstruction reports one: the run log-likelihood, the stop rule of optimize, the rate step and the α step
- **Combined posteriors**: the weight of category k at site s is P(k|s) = e^{ℓ_{s,k}} / Σ_j e^{ℓ_{s,j}}, and the posterior of a node state is Σ_k P(k|s) P_k(x|s), the marginal reconstruction under a rate mixture. The combined posteriors fill the node states of the existing types, so sequences, ambiguity codes, mutations, sampling and assembly read them without change. Leaf posteriors of ambiguous and missing states combine the category posteriors in the same way. A sparse combined node state keeps as explicit sites the union of the variable sites of the categories, plus the sites whose combined posterior differs from the combined group profile by more than the sparse profile tolerance (1e-4). Edge substitutions come from the combined node states
- **Edge contributions**: for edge e, the branch-length objective of the mixture is again a sum of exponentials in the branch length, with K * S terms: the coefficients of category k are the coefficients without gamma, built from the messages of pass k, with the eigenvalues r_k * λ_j, weighted by (1/K) e^{ℓ_{s,k} - log f_{s,k}(t_e) - m_s}, where f_{s,k}(t_e) is the site term at the branch length of the pass and m_s the largest exponent over k at site s. With K = 1 the contribution equals the one without gamma. The contribution and the full log-likelihood differ by a constant per edge, so tests compare differences
- **No zero-length certificate**: a mixture has several decay rates r_k * λ, also for JC69 and F81, so it gives no certified decay rate, and the branch-length optimizers do not accept a zero branch length from the derivative at zero alone

### Shape estimation and schedules

- **Rule**: α of each partition, or of each group that shares `gamma-shape`, maximizes the mixture log-likelihood (the sum over the group) with the bounded Brent solver of `argmin` in log α over [ln 0.02, ln 100] (`fn optimize_gamma_shape()`, [packages/treetime/src/gtr/gamma_shape.rs#L14](../../packages/treetime/src/gtr/gamma_shape.rs#L14)). Each evaluation rebuilds the category rates and runs the backward passes of all categories. The solver takes no start point; the first passes of a run use α = 1
- **Rate step**: the rate step of a mixture partition maximizes the mixture log-likelihood with the bracket and Brent search of the rate step without gamma (`fn optimize_mixture_rate()`, [packages/treetime/src/gtr/refinement.rs#L204](../../packages/treetime/src/gtr/refinement.rs#L204))
- **Optimize**: each iteration runs the branch-length step, the damping, the rate step, the α step, and the normalization (`fn update_partition_rates()`, [packages/treetime/src/optimize/partition_rates.rs#L13](../../packages/treetime/src/optimize/partition_rates.rs#L13)). A run with one partition has no rate step and runs the α step at the same place
- **Timetree**: each of the two ML branch-length passes before the dating rounds ends with the rate step, the α step and the normalization; the dating rounds keep α fixed
- **Ancestral**: one α step after the model fit, on the input branch lengths, before the final marginal update (`fn estimate_gamma_shapes()`, [packages/treetime/src/ancestral/plan.rs#L113](../../packages/treetime/src/ancestral/plan.rs#L113))

### Outputs

The model document (`--output-model`) gives each partition a `gamma` entry with the number of categories, α and the category rates; the entry is absent without `--gamma-categories`. The app shows α and K in the model table.

## Why

- **Rate variation across sites**: sites of one gene evolve at different speeds, and a model with one rate per partition underestimates long branches when slow and fast sites mix. Discrete gamma with four categories is the common approximation in IQ-TREE, RAxML-NG, BEAST and PhyML
- **One pass per category** keeps the propagation code of dense and sparse mode unchanged and keeps sparse fixed positions exact, because each pass uses one scalar rate. A fused pass with a category axis would need new message types in both representations
- **Mean-rate discretization** with equal weights matches the defaults of RAxML-NG and IQ-TREE, so the likelihoods of the three tools agree on the same model
- **Brent in log α**: the likelihood is smooth and usually unimodal in log α, the bracket is fixed by the shape range, and the solver needs no derivatives of the passes

## Verification

v0 has no gamma model, so the oracles are IQ-TREE, identities, and simulations with a known truth:

- **IQ-TREE**: on a simulated gap-free alignment with fixed branch lengths and model parameters, the per-site log-likelihoods of TreeTime equal those of IQ-TREE 3.1.4 (`GTR{...}+FU{...}+G4{α}`) for α of 0.2, 1 and 5: without gamma at 1e-12 in dense and sparse mode; with gamma at 1e-10 in dense mode, at 1e-8 in sparse mode for α of 1 and 5, and at 1e-6 in sparse mode for α of 0.2, whose slowest category rate (0.00053) meets the sparse profile tolerance. The category rates agree at 1e-8. IQ-TREE prints the per-site values of `-wsl` with six digits, so the capture runs IQ-TREE once per site with `-prec 17` (`packages/treetime/src/gtr/__tests__/test_gm_gamma_iqtree.rs`, fixture `gm_gamma_iqtree` in `dev/registry/reference-files.toml`)
- **Identities**: a single-category mixture equals the plain pass (dense profiles and contributions at 1e-12, sparse contributions bitwise); dense and sparse give the same mixture log-likelihood and posteriors; the difference of a mixture edge contribution between two branch lengths equals the difference of the full mixture log-likelihoods; the sparse per-site log-likelihoods sum to the pass log-likelihood; the category rates have mean 1 for every α in the range (`packages/treetime/src/partition/marginal/mixture/__tests__/`, `packages/treetime/src/gtr/__tests__/test_gamma.rs`)
- **Sparse tolerance**: the sparse state groups absorb resolved sites up to the profile tolerance of 1e-4, so sparse contribution differences deviate from full-pass differences by about 1e-5 on short branches, also with K = 1 and without gamma ([kb/issues/M-ancestral-dense-sparse-divergence.md](../issues/M-ancestral-dense-sparse-divergence.md)). The sparse exactness tests of the mixture use longer branches (lengths times 6) and α = 5
- **Shared shape**: a group of two partitions with identical data gives the α of one partition (`test_mixture_shared_shape.rs`)
