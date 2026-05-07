# Dense optimize iteration is slow

The `optimize.md` design document notes: "the iteration is slow for dense (which comes down to ancestral being compute heavy)."

Dense marginal reconstruction computes full probability profiles at every node for every alignment position. Each backward/forward pass costs $O(N \cdot s^2 \cdot L)$ where $N$ is the number of nodes, $s$ is the alphabet size, and $L$ is the alignment length -- the standard complexity of <a id="cite-1"></a>[Felsenstein 1981](https://doi.org/10.1007/BF01734359) [[1](#ref-1)] pruning. For a 30,000-site SARS-CoV-2 alignment with 500 taxa and $s = 5$, each pass processes $\sim 500 \times 25 \times 30{,}000 = 375{,}000{,}000$ operations. Sparse mode avoids this by operating only on variable sites (typically 1-5% of $L$).

## Current behavior

Each optimize iteration calls `update_marginal()` (backward + forward pass over the full dense profile matrix) followed by `run_optimize_mixed()` (per-edge Newton optimization). The marginal passes dominate wall time. The mpox dataset (`M-timetree-marginal-dense-mpox-slow`) is a known slow case for dense timetree inference.

## Potential optimizations

- **Column compression**: group identical alignment columns and weight by count, as in site pattern compression for likelihood evaluation. The marginal passes currently process every column individually. Compressed columns reduce $L$ to the number of unique patterns.
- **SIMD vectorization**: the inner loop multiplies profile vectors by transition probability matrices. AVX2/NEON intrinsics could process multiple sites or states per instruction, as demonstrated by the Phylogenetic Likelihood Library <a id="cite-2"></a>[Flouri et al. 2015](https://doi.org/10.1093/sysbio/syu084) [[2](#ref-2)].
- **Sparse fallback**: for alignments with low diversity, dense mode does redundant work on invariant columns. An adaptive mode could use dense for variable sites and a scalar shortcut for fixed sites.

## v0 comparison

v0 always uses dense representation. v0 mitigates the cost through column compression (grouping identical alignment columns by multiplicity). v1 dense does not compress columns, processing each individually. v1 sparse mode is the primary mitigation, reducing $L$ to variable sites only.

## Implementation pointers

Start with profiling (`cargo flamegraph` or `perf` via `./dev/docker/run`) on the mpox dataset in dense mode to confirm where time is spent. The expected hotspot is `marginal_backward`/`marginal_forward` in `marginal_dense.rs`.

Open questions for the implementer:

- Column compression (v0's approach) vs sparse-for-invariant-sites hybrid -- which gives better speedup for typical viral datasets? Column compression benefits high-redundancy alignments (SARS-CoV-2); the hybrid benefits low-diversity alignments.
- Is SIMD worth the maintenance cost? The inner loop is small ($s = 4$ or $5$ states). SIMD gains may be marginal compared to reducing $L$ via compression. PLL uses SIMD for $s = 20$ (amino acids) where the matrix multiply is larger.
- Could `ndarray` parallelism (via `rayon` or `ndarray-parallel`) parallelize across sites within a single pass? The per-site computations are independent.

The $O(N \cdot s^2 \cdot L)$ complexity is inherent -- the question is constant factors and whether $L$ can be reduced.

## Related

- [M-timetree-marginal-dense-mpox-slow](M-timetree-marginal-dense-mpox-slow.md) -- specific dataset where dense timetree is disproportionately slow
- [Command-relationships report](../reports/command-relationships/README.md)
- [Iterative tree refinement: status and recommendations](../reports/iterative-tree-refinement/10-status-and-recommendations.md)

## References

1. <a id="ref-1"></a> Felsenstein, Joseph. 1981. "Evolutionary Trees from DNA Sequences: A Maximum Likelihood Approach." _Journal of Molecular Evolution_ 17:368-376. https://doi.org/10.1007/BF01734359 [↩](#cite-1)
2. <a id="ref-2"></a> Flouri, Tomas, Fernando Izquierdo-Carrasco, Diego Darriba, Andre J. Aberer, Lam-Tung Nguyen, Bui Quang Minh, Arndt Von Haeseler, and Alexandros Stamatakis. 2015. "The Phylogenetic Likelihood Library." _Systematic Biology_ 64(2):356-362. https://doi.org/10.1093/sysbio/syu084 [↩](#cite-2)
