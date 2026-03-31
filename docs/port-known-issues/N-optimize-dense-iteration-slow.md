# Dense optimize iteration is slow

The `optimize.md` design document notes: "the iteration is slow for dense (which comes down to ancestral being compute heavy)."

Dense marginal reconstruction computes full probability profiles at every node for every alignment position. Each backward/forward pass costs $O(N \cdot s^2 \cdot L)$ where $N$ is the number of nodes, $s$ is the alphabet size, and $L$ is the alignment length. For a 30,000-site SARS-CoV-2 alignment with 500 taxa and $s = 5$, each pass processes $\sim 500 \times 25 \times 30{,}000 = 375{,}000{,}000$ operations. Sparse mode avoids this by operating only on variable sites (typically 1-5% of $L$).

## Current behavior

Each optimize iteration calls `update_marginal()` (backward + forward pass over the full dense profile matrix) followed by `run_optimize_mixed()` (per-edge Newton optimization). The marginal passes dominate wall time. The mpox dataset (`M-timetree-marginal-dense-mpox-slow`) is a known slow case for dense timetree inference.

## Potential optimizations

- **Column compression**: group identical alignment columns and weight by count, as in site pattern compression for likelihood evaluation. The marginal passes currently process every column individually. Compressed columns reduce $L$ to the number of unique patterns.
- **SIMD vectorization**: the inner loop multiplies profile vectors by transition probability matrices. AVX2/NEON intrinsics could process multiple sites or states per instruction.
- **Sparse fallback**: for alignments with low diversity, dense mode does redundant work on invariant columns. An adaptive mode could use dense for variable sites and a scalar shortcut for fixed sites.

## Scope

Profiling is needed to confirm where time is spent before implementing optimizations. The $O(N \cdot s^2 \cdot L)$ complexity is inherent to the Felsenstein pruning algorithm on dense profiles -- the question is whether constant factors can be reduced.

## Related

- [M-timetree-marginal-dense-mpox-slow](M-timetree-marginal-dense-mpox-slow.md) -- specific dataset where dense timetree is disproportionately slow
