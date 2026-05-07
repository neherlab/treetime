# Chapter 8: Initial branch length estimation

[Back to index](README.md) | Previous: [Chapter 7: Polytomy resolution](7-polytomy-resolution.md) | Next: [Chapter 9: The iteration loop](9-iteration-loop.md)

## Why initial guesses matter

The iteration loop ([Chapter 9](9-iteration-loop.md)) starts from initial branch lengths and refines them. The quality of the starting point affects convergence speed and correctness under multimodal likelihoods ([Chapter 5](5-branch-length-optimization.md)).

The input tree from the tree builder provides reasonable estimates. The initial guess step adjusts them using the alignment data before the first optimization iteration.

## The basic idea

The simplest initial guess is the fraction of positions where parent and child differ:

```
initial_bl = #differences / alignment_length
```

This is biased: multiple hits (two substitutions at the same position) appear as one difference. The true substitution count is higher than the observed difference count, with the bias growing for longer branches. For short branches (typical in viral phylogenetics), the bias is small.

## v0: Hamming distance with Brent bracket

v0 computes the Hamming distance between parent and child MAP sequences and divides by the alignment length:

```python
initial_bl = max(one_mutation, hamming_distance / alignment_length)
```

This seeds the bracket for Brent's method. Since Brent searches within `[0, 4.0]` regardless of starting point, the initial guess mainly affects convergence speed, not final accuracy.

v0 code: `optimal_t_compressed()` in [`packages/legacy/treetime/treetime/gtr.py#L816-L920`](../../../packages/legacy/treetime/treetime/gtr.py#L816-L920).

## v1: substitution count from marginal profiles

v1 counts substitutions from the marginal reconstruction and divides by the **effective alignment length** -- positions where both parent and child have canonical (non-gap, non-ambiguous) states:

```
initial_bl = #subs / effective_length
```

For sparse partitions, substitutions come from Fitch-reconstructed states. For dense partitions, substitutions come from MAP states (argmax of node posteriors). Gap positions are excluded from both numerator and denominator, giving correctly scaled rates for gappy alignments.

v1 code: `initial_guess_mixed()` in [`packages/treetime/src/commands/optimize/optimize_unified.rs#L729`](../../../packages/treetime/src/commands/optimize/optimize_unified.rs#L729).

## The GTR model interaction

The initial guess runs before the first optimization iteration but after the GTR model has been inferred. A subtlety: the marginal reconstruction (from which substitutions are counted) depends on the GTR model, and the GTR model is initially a dummy JC69. v1 runs marginal reconstruction twice -- once with dummy JC69, then once with the inferred GTR -- before computing the initial guess.

v1 code: the sequence in `run_optimize()` at [`packages/treetime/src/commands/optimize/run.rs`](../../../packages/treetime/src/commands/optimize/run.rs):

```
1. Create sparse/dense partitions with dummy JC69
2. Compress sequences (Fitch)
3. Infer GTR from compressed data
4. Replace dummy GTR with inferred GTR
5. Re-run marginal with inferred GTR
6. Compute initial_guess_mixed()      <-- uses fresh marginal profiles
7. Enter optimization loop
```

## References

- Felsenstein, J. 1981. "Evolutionary Trees from DNA Sequences." _J. Mol. Evol._ 17:368-376. https://doi.org/10.1007/BF01734359
- Sagulenko, P., V. Puller, and R. A. Neher. 2018. "TreeTime." _Virus Evol._ 4(1):vex042. https://doi.org/10.1093/ve/vex042
