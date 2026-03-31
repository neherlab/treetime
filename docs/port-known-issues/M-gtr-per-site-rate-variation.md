# Per-site rate variation not implemented

v1 uses a scalar substitution rate `mu` shared across all alignment positions. The design document (`docs/algorithms/sequence_evolution.md:87-89`) specifies per-site rate variation via a vector `mu^a`, where each site `a` evolves at its own rate while sharing the same eigendecomposition.

## Background

Among-site rate variation (ASRV) models the fact that different alignment positions evolve at different rates. Conserved positions (active sites, structural contacts) change slowly; variable positions (loops, surface residues) change fast. Without ASRV, branch length estimates are biased because the model assumes uniform rates across all sites.

The standard approach in phylogenetics is the discrete gamma model <a id="cite-1"></a>[Yang 1994](https://doi.org/10.1007/BF00160154) [[1](#ref-1)]: approximate the continuous gamma distribution of rates with K discrete categories (typically K=4), each with a rate multiplier. This is the "+Γ" suffix in model notation (e.g., "GTR+Γ4").

The design document describes a simpler approach: store a rate vector `mu^a` per site. The eigendecomposition is shared across all sites because only the rate scaling changes, not the rate matrix structure. The matrix exponential becomes `exp(Q * mu_a * t)` where `Q` is decomposed once and eigenvalues are scaled by `mu_a` per site.

## v1 current state

`packages/treetime/src/gtr/gtr.rs:173:` defines `pub mu: f64` as a scalar. All sites share this single rate. The `is_site_specific` field at `gtr.rs:169:` exists but is always `false`.

## v0 implementation

`GTR_site_specific` in `packages/legacy/treetime/treetime/gtr_site_specific.py` (495 lines) implements both per-site `mu` and per-site `pi`. The per-site `mu` path reuses the shared eigendecomposition, scaling eigenvalues per site.

## Implementation pointers

The core change is `mu: f64` to an enum or newtype that supports both scalar (all sites same rate) and per-site (vector) modes. A clean approach: keep `mu: f64` as the default scalar rate and add an optional `site_rates: Option<Array1<f64>>` field. When present, `expQt` multiplies eigenvalues by `site_rates[pos]` instead of scalar `mu`. The eigendecomposition stays shared.

v0's `GTR_site_specific` in `gtr_site_specific.py` stores per-site mu as a vector and overrides `expQt` to broadcast. Study that implementation for the scaling logic.

The main propagation concern is the marginal backward/forward passes: they call `expQt` per edge. With per-site rates, the transition matrix varies by site. Dense passes already iterate over sites (matrix-vector multiply per site), so the per-site rate slots in naturally. Sparse passes operate on variable sites only and would need the rate at each variable position.

The rate vector itself can be initialized from v0's approach: gamma-distributed rates with K=4 discrete categories, or empirical site rates from mutation counts. The estimation method is a separate concern from the `mu` representation change.

## Site-specific equilibrium frequencies

The design document (`docs/algorithms/sequence_evolution.md:85-86`) also describes a more general case: "If instead the equilibrium frequencies $\pi$ vary from site to site, then eigenvalues and eigenvectors change along the sequence." This means the matrix $e^{Q^a t}$ becomes site-specific, requiring per-site eigendecomposition rather than just per-site rate scaling.

Per-site $\mu$ (rate only): shared eigendecomposition, scaled eigenvalues. Computational cost: $O(n \cdot L \cdot s)$ where $s$ is alphabet size.

Per-site $\pi$ (equilibrium frequencies): per-site eigendecomposition. Computational cost: $O(n \cdot L \cdot s^2)$ or $O(n \cdot L \cdot s^3)$ depending on caching. The symmetrization trick ($\tilde{Q} = D^{-1} Q D$ where $D = \text{diag}(\sqrt{\pi})$) must be applied per site.

v1 has partial infrastructure for this: `gtr_site_specific.rs` implements per-site eigendecomposition but is not integrated into the partition system (tracked in [L-gtr-site-specific-partition-integration](L-gtr-site-specific-partition-integration.md)).

Context-dependent substitution models where both rates and equilibrium frequencies vary by position are described by <a id="cite-2"></a>[Siepel and Haussler 2004](https://doi.org/10.1093/molbev/msh039) [[2](#ref-2)].

These are two independent decisions:

- Per-site $\mu$: straightforward, reuses shared eigendecomposition (this issue)
- Per-site $\pi$: computationally expensive, per-site eigendecomposition required (separate scope)

## Related

- [L-gtr-site-specific-partition-integration](L-gtr-site-specific-partition-integration.md) - full site-specific GTR (per-site $\pi$) not yet integrated into partition system
- [docs/algorithms/sequence_evolution.md](../algorithms/sequence_evolution.md) - design document specifying per-site rate variation (lines 81-89)

## References

1. <a id="ref-1"></a> Yang, Ziheng. 1994. "Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites: Approximate Methods." _Journal of Molecular Evolution_ 39(3):306-314. https://doi.org/10.1007/BF00160154 [↩](#cite-1)
2. <a id="ref-2"></a> Siepel, Adam, and David Haussler. 2004. "Phylogenetic Estimation of Context-Dependent Substitution Rates by Maximum Likelihood." _Molecular Biology and Evolution_ 21(3):468-488. https://doi.org/10.1093/molbev/msh039 [↩](#cite-2)
