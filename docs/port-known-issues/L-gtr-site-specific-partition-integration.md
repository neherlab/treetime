# Site-specific GTR not integrated into partition system

The mathematical core for site-specific GTR models (`GTRSiteSpecific`) is implemented with full test coverage, but the production pipeline cannot use it because the partition system holds `pub gtr: GTR` (scalar model only).

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

- Change `pub gtr: GTR` to an enum or trait abstraction at 7 call sites across partition types (`PartitionMarginalDense`, `PartitionMarginalSparse`, `PartitionDiscrete`, `PartitionLikelihood`, `OptimizeDense`, and test input struct)
- The 2D `expQt()` callers need to handle the 3D site-specific case (different matrix per site)
- CLI flag `--site-specific-gtr` exists (hidden) on the ancestral command but returns an error
- Sequence compression (sparse representation) is incompatible with site-specific GTR because compressed positions share transition matrices but site-specific models give each position its own matrix

## Related

- Dead `is_site_specific: bool` field on `GTR` struct at `packages/treetime/src/gtr/gtr.rs:169:` should be removed when this integration is done
- Per-site rate variation (`M-gtr-per-site-rate-variation.md`) is a simpler feature where only mu varies per site (shared eigendecomposition). `GTRSiteSpecific` provides the full model (per-site mu AND per-site pi)
