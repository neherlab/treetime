# Stochastic polytomy event generator not implemented

v1 implements greedy deterministic polytomy resolution. v0 also supports a mutation-conditioned stochastic event generator for topology refinement.

## Background

A <a id="gloss-use-1"></a>polytomy <sup>[1](#gloss-1)</sup> is a node with more than two children. Resolving polytomies is optional topology refinement: timetree inference can retain a multifurcation, while resolution introduces additional internal branches and divergence times.

The v0 generator is not an ordinary Kingman sampler. In `generate_subtree()` at [`packages/legacy/treetime/treetime/treetime.py#L872-L1011`](../../packages/legacy/treetime/treetime/treetime.py#L872-L1011), it:

1. Sorts children by `time_before_present`.
2. Marks only mutation-free branches as ready to coalesce.
3. Samples competing mutation-removal and coalescence events using rates that depend on mutation counts and the ready subset.
4. Coalesces a sampled pair from the ready subset.
5. Stops when time is exhausted, potentially leaving a multifurcation.

The branch-merger rate can come from `merger_model`; mutation timing uses `gtr.mu * L`. Random draws use NumPy's `default_rng` through `self.rng.exponential`, `self.rng.random`, and `self.rng.choice` in event-dependent call order.

## v1 status

v1 uses greedy pairwise merging with Brent optimization in [`packages/treetime/src/timetree/optimization/polytomy.rs`](../../packages/treetime/src/timetree/optimization/polytomy.rs). It has no counterpart for v0's stochastic event generator.

The public tracker discusses information sources and efficiency for polytomy resolution [[issue](https://github.com/neherlab/treetime/issues/109)] and comparison of stochastic and greedy resolution [[issue](https://github.com/neherlab/treetime/issues/313)]. This issue tracks the narrower v0 parity gap.

## Impact

- v1 cannot reproduce v0's stochastic topology-generation behavior.
- Greedy and stochastic resolution can produce different topology and timing inputs for downstream inference.
- v0's specialized event generator is not a calibrated posterior over resolutions. The star tree paradox concerns posterior support for arbitrary resolved topologies when data arise from a star tree; proposed remedies give unresolved topologies explicit prior mass rather than randomly resolving them <a id="cite-1"></a>[Lewis et al. 2005](https://doi.org/10.1080/10635150590924208) [[1](#ref-1)]. A calibrated topology-uncertainty feature would require a separately approved target distribution and validation contract.

## Algorithm detail

See [kb/algo/unimplemented.md](../algo/unimplemented.md#stochastic-polytomy-resolution) for the v0 algorithm walkthrough.

## Decisions required

- Choose exact seeded parity, which requires the pinned NumPy reference version, NumPy `default_rng` PCG64, parity with NumPy's exponential, uniform-random, and choice transformations, and identical draw ordering; or choose a Rust RNG with distributional and invariant-based comparison.
- Define whether the CLI exposes v0's specialized generator as parity behavior or adopts a separately approved stochastic topology model.

No implementation ticket is ready until these contracts are approved.

## Cross-references

- [kb/reports/iterative-tree-refinement/7-polytomy-resolution.md](../reports/iterative-tree-refinement/7-polytomy-resolution.md)

## Glossary

1. <a id="gloss-1"></a> **Polytomy.** A tree node with more than two children. [↩](#gloss-use-1)

## References

1. <a id="ref-1"></a> Lewis, Paul O., Mark T. Holder, and Kent E. Holsinger. 2005. "Polytomies and Bayesian phylogenetic inference." _Systematic Biology_ 54(2):241-253. https://doi.org/10.1080/10635150590924208 [↩](#cite-1)
