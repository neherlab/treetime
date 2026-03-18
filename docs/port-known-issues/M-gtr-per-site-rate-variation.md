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

## References

- Yang Z (1994). "Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods." J Mol Evol 39:306-314.
- `docs/algorithms/sequence_evolution.md` lines 87-89.
