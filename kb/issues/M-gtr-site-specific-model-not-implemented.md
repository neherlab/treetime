# Site-specific GTR is not implemented

v1 has no site-specific substitution model. The hidden `--site-specific-gtr` flag of `ancestral` and `homoplasy` ([packages/app-commands/src/commands/ancestral/args.rs#L305-L311](../../packages/app-commands/src/commands/ancestral/args.rs#L305-L311), [packages/app-commands/src/commands/homoplasy/args.rs#L169-L175](../../packages/app-commands/src/commands/homoplasy/args.rs#L169-L175)) is accepted and then rejected with the error `--site-specific-gtr is not implemented` [`packages/treetime/src/ancestral/plan.rs#L30-L35`](../../packages/treetime/src/ancestral/plan.rs#L30-L35). v0 supports the model in production, so this is a parity gap.

## v0 behavior

`treeanc.infer_gtr(site_specific=True)` at `packages/legacy/treetime/treetime/treeanc.py:1500-1632` creates a `GTR_site_specific` model and uses it transparently through the same tree traversal code. `GTR_site_specific` in `packages/legacy/treetime/treetime/gtr_site_specific.py` holds per-site equilibrium frequencies $\pi^a$ and rates $\mu^a$ and infers them iteratively from per-site mutation counts.

## Work needed

- A site-specific model with per-site eigendecomposition, transition matrices `[n_states, n_states, seq_len]`, and backward and forward propagation
- Inference of $W$, $\pi^a$ and $\mu^a$ from per-site mutation counts, with golden-master tests against v0
- Model access that both scalar and site-specific models can honor: partitions and reconstructions store the concrete scalar `GTR` type, for example [`packages/treetime/src/partition/create.rs#L95`](../../packages/treetime/src/partition/create.rs#L95) and [`packages/treetime/src/partition/optimize/dense.rs#L17`](../../packages/treetime/src/partition/optimize/dense.rs#L17)
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
> **Decision required.** With a site-specific model and sparse mode requested, the options are: reject the combination with an error, or select dense mode automatically. The flag help already states that the model requires `--dense true` and is incompatible with the sparse representation.

Study v0's `treeanc.infer_gtr(site_specific=True)` path for the inference-to-traversal integration pattern. v0 uses duck typing (both `GTR` and `GTR_site_specific` expose `expQt`), which maps naturally to a Rust trait.

## Related

- [M-gtr-per-site-rate-variation.md](M-gtr-per-site-rate-variation.md) - simpler feature where only $\mu$ varies per site (shared eigendecomposition)
- [kb/algo/unimplemented.md](../algo/unimplemented.md) - site-specific GTR in v0

## References

1. <a id="ref-1"></a> Yang, Ziheng. 1994. "Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites: Approximate Methods." _Journal of Molecular Evolution_ 39(3):306-314. https://doi.org/10.1007/BF00160154 [↩](#cite-1)
2. <a id="ref-2"></a> Siepel, Adam, and David Haussler. 2004. "Phylogenetic Estimation of Context-Dependent Substitution Rates by Maximum Likelihood." _Molecular Biology and Evolution_ 21(3):468-488. https://doi.org/10.1093/molbev/msh039 [↩](#cite-2)
