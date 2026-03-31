# Grid search fallback covers narrow range at zero branch length

The grid search fallback at [packages/treetime/src/commands/optimize/optimize_unified.rs#L275](../../packages/treetime/src/commands/optimize/optimize_unified.rs#L275) uses a range proportional to the current branch length:

```rust
let branch_lengths = ndarray::Array1::linspace(0.1 * one_mutation, 1.5 * branch_length + one_mutation, 100);
```

When `branch_length = 0.0`, this becomes `linspace(0.1/L, 1.0/L, 100)` where L is total alignment length. For a 1000-site alignment, the grid covers only `[1e-4, 1e-3]`, missing branches longer than `1/L` substitutions per site.

## Scientific background

The phylogenetic likelihood surface for branch lengths can be multimodal under general substitution models <a id="cite-1"></a>[Dinh and Matsen 2017](https://doi.org/10.1214/16-AAP1240) [[1](#ref-1)]. The grid search exists as a fallback when Newton-Raphson fails on non-concave surfaces. For the grid to serve its purpose, it must cover the biologically plausible range of branch lengths. Branch lengths in phylogenetic trees span 3-4 orders of magnitude ($10^{-5}$ to $10^{-1}$ subs/site), so a grid restricted to one decade near $1/L$ misses the long-branch regime entirely.

## Impact

Medium. The outer coordinate-ascent loop partially compensates by re-optimizing each branch in subsequent iterations, but convergence is slower. Edges with zero initial length that should converge to long branches (e.g. after tree rearrangement) are trapped near zero until the slow re-expansion of the grid range catches up.

## v0 comparison

v0 uses `scipy.optimize.minimize_scalar(method='bounded')` (Brent's method) with fixed bracket $[0, 4.0]$ in $\sqrt{t}$ space <a id="cite-2"></a>[Sagulenko, Puller, and Neher 2018](https://doi.org/10.1093/ve/vex042) [[2](#ref-2)]. The bracket does not depend on the current branch length, avoiding this issue. Brent's method is derivative-free and does not require a grid <a id="cite-3"></a>[Brent 1973](https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf) [[3](#ref-3)].

## Cross-links

- [M-optimize-grid-zero-ignores-indels.md](M-optimize-grid-zero-ignores-indels.md) -- related grid search issue affecting indel-aware zero comparison
- [H-optimize-sparse-hessian-multiplicity.md](H-optimize-sparse-hessian-multiplicity.md) -- Hessian bug causes more edges to fall through to grid search
- [docs/reports/optimization-methods/2-branch-length-optimization.md](../../docs/reports/optimization-methods/2-branch-length-optimization.md) -- comparison table of fallback strategies across tools

## Proposed fix

Use a minimum upper bound independent of the current branch length:

```rust
let upper = f64::max(1.5 * branch_length + one_mutation, 0.5);
let branch_lengths = ndarray::Array1::linspace(0.1 * one_mutation, upper, 100);
```

A log-spaced grid is a stronger alternative. Branch lengths span orders of magnitude, and a linear grid wastes resolution on the narrow range near zero while under-sampling the long-branch tail. A log-spaced grid covers the full biologically plausible range with uniform resolution per decade:

```rust
let log_min = (0.1 * one_mutation).ln();
let log_max = 0.5_f64.ln();
let branch_lengths = ndarray::Array1::linspace(log_min, log_max, 100).mapv(f64::exp);
```

For reference, RAxML uses derivative-based optimization with Brent fallback on bounded intervals that do not shrink with the current estimate <a id="cite-4"></a>[Stamatakis 2014](https://doi.org/10.1093/bioinformatics/btu033) [[4](#ref-4)].

## References

1. <a id="ref-1"></a> Dinh, Vu, and Frederick A. Matsen IV. 2017. "The Shape of the One-Dimensional Phylogenetic Likelihood Function." _Annals of Applied Probability_ 27(3):1432-1461. https://doi.org/10.1214/16-AAP1240 [↩](#cite-1)
2. <a id="ref-2"></a> Sagulenko, Pavel, Vadim Puller, and Richard A. Neher. 2018. "TreeTime: Maximum-Likelihood Phylodynamic Analysis." _Virus Evolution_ 4(1):vex042. https://doi.org/10.1093/ve/vex042 [↩](#cite-2)
3. <a id="ref-3"></a> Brent, Richard P. 1973. "Algorithms for Minimization without Derivatives." Prentice-Hall. https://maths-people.anu.edu.au/~brent/pd/rpb011i.pdf [↩](#cite-3)
4. <a id="ref-4"></a> Stamatakis, Alexandros. 2014. "RAxML Version 8: A Tool for Phylogenetic Analysis and Post-Analysis of Large Phylogenies." _Bioinformatics_ 30(9):1312-1313. https://doi.org/10.1093/bioinformatics/btu033 [↩](#cite-4)
