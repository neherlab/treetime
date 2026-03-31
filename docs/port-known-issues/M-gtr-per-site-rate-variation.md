# Per-site rate variation not implemented

v1 uses a scalar substitution rate `mu` shared across all alignment positions. The design document (`docs/algorithms/sequence_evolution.md:87-89`) specifies per-site rate variation via a vector `mu^a`, where each site `a` evolves at its own rate while sharing the same eigendecomposition.

## Background

Among-site rate variation (ASRV) models the fact that different alignment positions evolve at different rates. Conserved positions (active sites, structural contacts) change slowly; variable positions (loops, surface residues) change fast. Without ASRV, branch length estimates are biased because the model assumes uniform rates across all sites.

The standard approach in phylogenetics is the discrete gamma model (Yang 1994): approximate the continuous gamma distribution of rates with K discrete categories (typically K=4), each with a rate multiplier. This is the "+Γ" suffix in model notation (e.g., "GTR+Γ4").

The design document describes a simpler approach: store a rate vector `mu^a` per site. The eigendecomposition is shared across all sites because only the rate scaling changes, not the rate matrix structure. The matrix exponential becomes `exp(Q * mu_a * t)` where `Q` is decomposed once and eigenvalues are scaled by `mu_a` per site.

## v1 current state

`packages/treetime/src/gtr/gtr.rs:173:` defines `pub mu: f64` as a scalar. All sites share this single rate. The `is_site_specific` field at `gtr.rs:169:` exists but is always `false`.

## v0 implementation

`GTR_site_specific` in `packages/legacy/treetime/treetime/gtr_site_specific.py` (495 lines) implements both per-site `mu` and per-site `pi`. The per-site `mu` path reuses the shared eigendecomposition, scaling eigenvalues per site.

## Implementation path

Change `mu: f64` to `mu: Array1<f64>` (or keep scalar as default, vector as opt-in). In `expQt`, scale eigenvalues by `mu[site]` instead of scalar `mu`. The eigendecomposition (`eigenvalues`, `v`, `v_inv`) remains shared. This affects `gtr.rs` and all call sites that use `mu` for rate scaling.

## Site-specific equilibrium frequencies

The design document (`docs/algorithms/sequence_evolution.md:85-86`) also describes a more general case: "If instead the equilibrium frequencies $\pi$ vary from site to site, then eigenvalues and eigenvectors change along the sequence." This means the matrix $e^{Q^a t}$ becomes site-specific, requiring per-site eigendecomposition rather than just per-site rate scaling.

Per-site $\mu$ (rate only): shared eigendecomposition, scaled eigenvalues. Computational cost: $O(n \cdot L \cdot s)$ where $s$ is alphabet size.

Per-site $\pi$ (equilibrium frequencies): per-site eigendecomposition. Computational cost: $O(n \cdot L \cdot s^2)$ or $O(n \cdot L \cdot s^3)$ depending on caching. The symmetrization trick ($\tilde{Q} = D^{-1} Q D$ where $D = \text{diag}(\sqrt{\pi})$) must be applied per site.

v1 has partial infrastructure for this: `gtr_site_specific.rs` implements per-site eigendecomposition but is not integrated into the partition system (tracked in `L-gtr-site-specific-partition-integration`).

These are two independent decisions:

- Per-site $\mu$: straightforward, reuses shared eigendecomposition (this issue)
- Per-site $\pi$: computationally expensive, per-site eigendecomposition required (separate scope)

## References

- Yang Z (1994). "Maximum Likelihood Phylogenetic Estimation from DNA Sequences with Variable Rates over Sites: Approximate Methods." J Mol Evol 39(3):306-314. https://doi.org/10.1007/BF00160154
- Siepel A, Haussler D (2004). "Phylogenetic Estimation of Context-Dependent Substitution Rates by Maximum Likelihood." Mol Biol Evol 21(3):468-488. https://doi.org/10.1093/molbev/msh039 -- context-dependent models where both rates and equilibrium frequencies vary by position.
- `docs/algorithms/sequence_evolution.md` lines 81-86.
