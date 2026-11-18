# Site-specific GTR is not implemented

v1 has no site-specific substitution model, and no command has a flag for it (v0: `--site-specific-gtr` of `ancestral` and `homoplasy`) until the model exists. v0 supports the model in production, so this is a parity gap.

## v0 behavior

`treeanc.infer_gtr(site_specific=True)` at `packages/legacy/treetime/treetime/treeanc.py:1500-1632` creates a `GTR_site_specific` model and uses it transparently through the same tree traversal code. `GTR_site_specific` in `packages/legacy/treetime/treetime/gtr_site_specific.py` holds per-site equilibrium frequencies $\pi^a$ and rates $\mu^a$ and infers them iteratively from per-site mutation counts.

## Work needed

- A site-specific model with per-site eigendecomposition, transition matrices `[n_states, n_states, seq_len]`, and backward and forward propagation
- The `--site-specific-gtr` flag of `ancestral` and `homoplasy`, added together with the model
- Inference of $W$, $\pi^a$ and $\mu^a$ from per-site mutation counts, with golden-master tests against v0
- Model access that both scalar and site-specific models can honor: partitions and reconstructions store the concrete scalar `GTR` type, for example [`packages/treetime/src/partition/create.rs#L121`](../../packages/treetime/src/partition/create.rs#L121) and [`packages/treetime/src/partition/optimize/dense.rs#L8`](../../packages/treetime/src/partition/optimize/dense.rs#L8)
- Callers of the 2D `expQt()` that handle a different matrix per site
- A policy for sparse mode: compressed positions share transition matrices, while a site-specific model gives each position its own matrix

## Scientific background

Site-specific substitution models generalize the standard GTR framework by allowing model parameters to vary across alignment positions. The simplest form is among-site rate variation, where a per-site rate multiplier $\mu^a$ scales a shared rate matrix <a id="cite-1"></a>[Yang 1994](https://doi.org/10.1007/BF00160154) [[1](#ref-1)]. The full site-specific model allows both rates and equilibrium frequencies $\pi^a$ to vary per site, requiring per-site eigendecomposition of the rate matrix. <a id="cite-2"></a>[Siepel and Haussler 2004](https://doi.org/10.1093/molbev/msh039) [[2](#ref-2)] describe context-dependent substitution models with per-site eigendecomposition, showing improved fit for protein-coding regions where different codon positions have distinct substitution patterns.

For each site $a$, the rate matrix $Q^a$ is built from site-specific $\pi^a$ and $\mu^a$, symmetrized via $\tilde{Q}^a = D_a^{-1} Q^a D_a$ where $D_a = \text{diag}(\sqrt{\pi^a})$, and eigendecomposed independently.

## Design axes

### Model abstraction

- O1. Store an enum covering scalar and site-specific models. This makes supported variants exhaustive but forces consumers to handle dimensional differences at the enum boundary.
- O2. Define role-specific capabilities for transition propagation, likelihood evaluation, and rate refinement. This prevents consumers from depending on model operations they do not use, but associated output dimensionality must remain explicit.

> [!IMPORTANT]
> **Decision required.** The concrete `GTR` fields cannot hold a site-specific model, and changing them requires an approved architecture decision. The options are O1 (an enum of scalar and site-specific models) and O2 (role-specific capabilities). Evidence: partitions, reconstructions, and the branch length objective store the concrete `GTR` type (see Work needed), and v0 relies on duck typing between `GTR` and `GTR_site_specific` (see below).

The sparse incompatibility is a hard constraint: a site-specific model requires dense position-specific propagation.

> [!IMPORTANT]
> **Decision required.** With a site-specific model and sparse mode requested, the options are: reject the combination with an error, or select dense mode automatically.

Study v0's `treeanc.infer_gtr(site_specific=True)` path for the inference-to-traversal integration pattern. v0 uses duck typing (both `GTR` and `GTR_site_specific` expose `expQt`), which maps naturally to a Rust trait.

## Related

- [kb/decisions/discrete-gamma-rate-variation.md](../decisions/discrete-gamma-rate-variation.md) - rate variation across sites with one shared rate matrix (discrete gamma)
- [kb/algo/unimplemented.md](../algo/unimplemented.md) - site-specific GTR in v0

## References

1. <a id="ref-1"></a> Yang, Ziheng. 1994. "Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites: Approximate Methods." _Journal of Molecular Evolution_ 39(3):306-314. https://doi.org/10.1007/BF00160154 [↩](#cite-1)
2. <a id="ref-2"></a> Siepel, Adam, and David Haussler. 2004. "Phylogenetic Estimation of Context-Dependent Substitution Rates by Maximum Likelihood." _Molecular Biology and Evolution_ 21(3):468-488. https://doi.org/10.1093/molbev/msh039 [↩](#cite-2)

## Proposed implementation plan

No Rust source mentions a site-specific model or flag; [`packages/treetime/src/partition/create.rs#L121`](../../packages/treetime/src/partition/create.rs#L121) returns `Vec<GTR>` and [`packages/treetime/src/partition/optimize/dense.rs#L8`](../../packages/treetime/src/partition/optimize/dense.rs#L8) takes `&GTR`. The plan stops at the two open decisions above.

1. Get the user's choice of model abstraction (O1 or O2) and of the sparse policy, and record them in a new `kb/decisions/` file with approval.
2. Model: add `GtrSiteSpecific` in `packages/treetime/src/gtr/` with per-site `pi` (`Array2<f64>` `[n_states, seq_len]`), `mu` (`Array1<f64>`), shared `W`, per-site eigendecomposition, and `exp_qt` returning `Array3<f64>` `[n_states, n_states, seq_len]`, following `gtr_site_specific.py`.
   - O1: an enum `SubstitutionModel { Scalar(GTR), SiteSpecific(GtrSiteSpecific) }` replaces `GTR` in partitions, reconstructions, and the branch length objective.
   - O2: traits for propagation, likelihood, and rate refinement, implemented by both models; the stored type becomes generic or boxed per role.
3. Inference: port `GTR_site_specific.infer` (iterative per-site `pi`, `mu`) from per-site counts that `accumulate_mutation_counts` extends with a site axis.
4. Propagation: update dense marginal backward and forward passes to apply a per-site matrix; apply the sparse policy chosen in step 1.
5. CLI: add `--site-specific-gtr` to `ancestral` and `homoplasy`, regenerate schemas (`just gen`), and remove the 0.x message for it from `packages/app-cli/src/cli/v0_args.rs` if it exists.

Tests: golden master `test_gm_infer_gtr_site_specific` against v0 `treeanc.infer_gtr(site_specific=True)` (capture script imports v0) on `data/flu/h3n2/20` and `data/ebola/20` at 1e-6; unit test that a site-specific model with equal per-site `pi` and `mu` gives the same likelihood as the scalar `GTR` (oracle: the scalar model); dense marginal ancestral output vs v0 `ancestral --site-specific-gtr`.

Validation: `just lint-rs`, `just test-rs`, `just generated-check`, `just smoke` on the affected commands.
