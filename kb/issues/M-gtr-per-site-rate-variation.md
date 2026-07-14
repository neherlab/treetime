# Per-site rate variation is not inferred or wired end to end

v1 represents and propagates an optional per-site rate vector, but no command constructs that vector from input or inference and optimizer coverage is incomplete. The design document specifies a vector $\mu^a$, where $a$ indexes sites and each site evolves at its own supplied rate. [`kb/_raw/sequence_evolution.md#L87-L89`](../_raw/sequence_evolution.md#L87-L89)

## Background

<a id="gloss-use-1"></a>Among-site rate variation (ASRV) <sup>[1](#gloss-1)</sup> models the fact that different alignment positions evolve at different rates. Conserved positions (active sites, structural contacts) change slowly; variable positions (loops, surface residues) change fast. Without ASRV, branch length estimates are biased because the model assumes uniform rates across all sites.

The standard approach in phylogenetics is the discrete gamma model <a id="cite-1"></a>[Yang 1994](https://doi.org/10.1007/BF00160154) [[1](#ref-1)]: approximate the continuous gamma distribution of rates with K discrete categories (typically K=4), each with a rate multiplier. This is the "+Γ" suffix in model notation (e.g., "GTR+Γ4").

The design document's fixed rate-vector contract is distinct from discrete-gamma ASRV. For site $a$, a discrete-gamma model evaluates likelihood $L_a=\sum_k w_kL_a(r_k)$ over rate categories $k$, with category weights $w_k$ and rate multipliers $r_k$; a supplied vector evaluates one rate $\mu_a$ per site. Full site-specific GTR is a third contract in which equilibrium frequencies $\pi^a$ or other rate-matrix parameters also vary by site.

## v1 current state

`struct GTR` stores scalar `mu` together with optional `site_rates`. [`packages/treetime/src/gtr/gtr.rs#L184-L193`](../../packages/treetime/src/gtr/gtr.rs#L184-L193) Its transition-matrix paths apply the site-rate vector, and dense marginal propagation passes position-specific matrices. [`packages/treetime/src/gtr/gtr.rs#L303-L420`](../../packages/treetime/src/gtr/gtr.rs#L303-L420) [`packages/treetime/src/partition/marginal_passes.rs#L149-L158`](../../packages/treetime/src/partition/marginal_passes.rs#L149-L158)

The remaining gaps are rate-vector construction or loading, command wiring, compressed fixed-position behavior, optimizer coverage, validation, serialization, and end-to-end comparison with v0. Sparse explicit positions already index `site_rates[pos]`; this does not define how compressed fixed positions contribute to likelihood or optimization.

## v0 implementation

`GTR_site_specific` in `packages/legacy/treetime/treetime/gtr_site_specific.py` implements both per-site $\mu$ and per-site $\pi$. The reference loops over sites and recomputes the transition-matrix decomposition even when $\pi$ is shared. V1 can reuse a mathematically shared decomposition only when doing so preserves the approved output-parity contract.

## Decisions required

- Choose whether the target is a supplied fixed rate vector, discrete-gamma latent categories, v0 full site-specific inference, or separately exposed contracts, and define normalization for each approved contract.
- Specify how compressed fixed sites contribute to propagation and optimization; explicit sparse positions already map by alignment position.
- Define serialization, validation, and CLI behavior.
- Establish dense, sparse, and optimizer golden masters against v0, with independent analytical tests for any implementation strategy that differs internally.

No implementation ticket is ready until these contracts are decided.

## Site-specific equilibrium frequencies

The design document (`../_raw/sequence_evolution.md:85-86`) also describes a more general case: "If instead the equilibrium frequencies $\pi$ vary from site to site, then eigenvalues and eigenvectors change along the sequence." This means the matrix $e^{Q^a t}$ becomes site-specific, requiring per-site eigendecomposition rather than just per-site rate scaling.

Per-site $\mu$ (rate only): shared eigendecomposition, scaled eigenvalues. Computational cost: $O(n \cdot L \cdot s)$ where $s$ is alphabet size.

Per-site $\pi$ (equilibrium frequencies): per-site eigendecomposition. Computational cost: $O(n \cdot L \cdot s^2)$ or $O(n \cdot L \cdot s^3)$ depending on caching. The symmetrization trick ($\tilde{Q} = D^{-1} Q D$ where $D = \text{diag}(\sqrt{\pi})$) must be applied per site.

v1 has partial infrastructure for this: `gtr_site_specific.rs` implements per-site eigendecomposition but is not integrated into the partition system (tracked in [L-gtr-site-specific-partition-integration](N-gtr-site-specific-partition-integration.md)).

Context-dependent substitution models where both rates and equilibrium frequencies vary by position are described by <a id="cite-2"></a>[Siepel and Haussler 2004](https://doi.org/10.1093/molbev/msh039) [[2](#ref-2)].

These are two independent decisions:

- Per-site $\mu$: straightforward, reuses shared eigendecomposition (this issue)
- Per-site $\pi$: computationally expensive, per-site eigendecomposition required (separate scope)

## Related

- [L-gtr-site-specific-partition-integration](N-gtr-site-specific-partition-integration.md) - full site-specific GTR (per-site $\pi$) not yet integrated into partition system
- [../_raw/sequence_evolution.md](../_raw/sequence_evolution.md) - design document specifying per-site rate variation (lines 81-89)

## Glossary

1. <a id="gloss-1"></a> **Among-site rate variation (ASRV).** Variation in the evolutionary rate among alignment positions. [↩](#gloss-use-1)

## References

1. <a id="ref-1"></a> Yang, Ziheng. 1994. "Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: Approximate methods." _Journal of Molecular Evolution_ 39(3):306-314. https://doi.org/10.1007/BF00160154 [↩](#cite-1)
2. <a id="ref-2"></a> Siepel, Adam, and David Haussler. 2004. "Phylogenetic estimation of context-dependent substitution rates by maximum likelihood." _Molecular Biology and Evolution_ 21(3):468-488. https://doi.org/10.1093/molbev/msh039 [↩](#cite-2)
