# Site-specific GTR not integrated into partition system

The mathematical core for site-specific GTR models (`GTRSiteSpecific`) is implemented, but the production pipeline cannot use it because the partition capability `HasGtr` returns the concrete scalar `GTR` type [`packages/treetime/src/partition/traits.rs#L27`](../../packages/treetime/src/partition/traits.rs#L27). The abstraction names a capability while fixing one representation.

## Current state

`GTRSiteSpecific` at `packages/treetime/src/gtr/gtr_site_specific.rs` provides:

- Per-site eigendecomposition with per-site pi and mu
- `expQt()` returning 3D transition matrices `[n_states, n_states, seq_len]`
- `propagate_profile()` and `evolve()` for message passing
- Iterative inference from per-site mutation counts
- Linear interpolation for fast evaluation
- Golden-master validation against v0 at 1e-10

v0's production path: `treeanc.infer_gtr(site_specific=True)` at `packages/legacy/treetime/treetime/treeanc.py:1500-1632:` creates a `GTR_site_specific` model and uses it transparently through the same tree traversal code.

## Work needed

- Replace concrete scalar-model access at consumers with model capabilities that both scalar and site-specific implementations can honor.
- The 2D `expQt()` callers need to handle the 3D site-specific case (different matrix per site)
- CLI flag `--site-specific-gtr` exists (hidden) on the ancestral command but returns an error
- Sequence compression (sparse representation) is incompatible with site-specific GTR because compressed positions share transition matrices but site-specific models give each position its own matrix

## Scientific background

Site-specific substitution models generalize the standard GTR framework by allowing model parameters to vary across alignment positions. The simplest form is among-site rate variation, where a per-site rate multiplier $\mu^a$ scales a shared rate matrix <a id="cite-1"></a>[Yang 1994](https://doi.org/10.1007/BF00160154) [[1](#ref-1)]. The full site-specific model allows both rates and equilibrium frequencies $\pi^a$ to vary per site, requiring per-site eigendecomposition of the rate matrix. <a id="cite-2"></a>[Siepel and Haussler 2004](https://doi.org/10.1093/molbev/msh039) [[2](#ref-2)] describe context-dependent substitution models with per-site eigendecomposition, showing improved fit for protein-coding regions where different codon positions have distinct substitution patterns.

The `GTRSiteSpecific` implementation in v1 follows the per-site eigendecomposition approach: for each site $a$, the rate matrix $Q^a$ is constructed from site-specific $\pi^a$ and $\mu^a$, symmetrized via $\tilde{Q}^a = D_a^{-1} Q^a D_a$ where $D_a = \text{diag}(\sqrt{\pi^a})$, and eigendecomposed independently.

## Design axes

### Model abstraction

- O1. Store an enum covering scalar and site-specific models. This makes supported variants exhaustive but forces consumers to handle dimensional differences at the enum boundary.
- O2. Define role-specific capabilities for transition propagation, likelihood evaluation, and rate refinement. This prevents consumers from depending on model operations they do not use, but associated output dimensionality must remain explicit.

No option is selected. The current `HasGtr` contract cannot represent the implemented site-specific model, and changing that contract requires an approved architecture decision.

The sparse incompatibility (compressed positions share transition matrices) is a hard constraint: site-specific GTR requires dense position-specific propagation. Whether the command rejects sparse mode or selects dense mode automatically is a separate user-facing policy decision.

Study v0's `treeanc.infer_gtr(site_specific=True)` path for the inference-to-traversal integration pattern. v0 uses duck typing (both `GTR` and `GTR_site_specific` expose `expQt`), which maps naturally to a Rust trait.

## Related

- The `is_site_specific: bool` field formerly on `GTR` has been removed
- [M-gtr-per-site-rate-variation](M-gtr-per-site-rate-variation.md) - simpler feature where only $\mu$ varies per site (shared eigendecomposition); `GTRSiteSpecific` provides the full model (per-site $\mu$ AND per-site $\pi$)
- [../_raw/sequence_evolution.md](../_raw/sequence_evolution.md) - design document specifying site-specific models

## References

1. <a id="ref-1"></a> Yang, Ziheng. 1994. "Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites: Approximate Methods." _Journal of Molecular Evolution_ 39(3):306-314. https://doi.org/10.1007/BF00160154 [↩](#cite-1)
2. <a id="ref-2"></a> Siepel, Adam, and David Haussler. 2004. "Phylogenetic Estimation of Context-Dependent Substitution Rates by Maximum Likelihood." _Molecular Biology and Evolution_ 21(3):468-488. https://doi.org/10.1093/molbev/msh039 [↩](#cite-2)
