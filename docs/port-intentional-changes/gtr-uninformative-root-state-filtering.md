# Uninformative root_state filtering in GTR inference

This document describes an intentional deviation from v0 in the dense GTR inference path. v1 filters out alignment positions with uninformative marginal profiles when computing `root_state`, while v0 includes all positions regardless of signal content.

The change affects `root_state` computation in `get_mutation_counts_dense()` (`#get_mutation_counts_dense`) at `packages/treetime/src/gtr/infer_gtr/dense.rs:139-154:`. v0's equivalent code is at `packages/legacy/treetime/treetime/treeanc.py:1608-1613:`. Downstream, this affects equilibrium frequencies (pi), exchangeability matrix (W), and rate scalar (mu) in the GTR model.

Datasets with gap-only columns are affected. Measured impact: lassa_L_50 showed ~29% shift in W matrix elements before the golden master capture script was updated to match v1's filtering. After updating the capture script, all seven real datasets pass at 1e-6 tolerance.

## Background: GTR model and equilibrium frequencies

The General Time Reversible (GTR) model, introduced by Tavaré (1986), describes nucleotide substitution as a continuous-time Markov chain. It is the most general neutral, independent, finite-sites, time-reversible model for DNA sequence evolution. The model has two core components:

- **Equilibrium frequencies (pi)**: A vector where `pi[k]` is the stationary probability of state k. At equilibrium, the fraction of positions in state k converges to `pi[k]`. For DNA, this is a 4-element vector (pi_A, pi_G, pi_C, pi_T) summing to 1.
- **Exchangeability matrix (W)**: A symmetric matrix where `W[i,j]` describes the relative rate of exchange between states i and j. For DNA, this has 6 independent parameters representing rates between each nucleotide pair.

The instantaneous rate matrix Q combines these: `Q[i,j] = W[i,j] * pi[j]` for i != j, with diagonal elements set so rows sum to zero. The time-reversibility condition requires `pi[i] * Q[i,j] = pi[j] * Q[j,i]` (detailed balance), which allows likelihood calculations without knowing the root position (Tavaré, 1986).

Phylogenetic inference estimates these parameters from observed substitution patterns. The standard approach:

1. Count observed substitutions between states (nij matrix)
2. Estimate time spent in each state (Ti vector)
3. Use the root sequence composition as a prior on pi

The iterative solver in both v0 (`packages/legacy/treetime/treetime/gtr.py:492-599:`) and v1 (`packages/treetime/src/gtr/infer_gtr/common.rs:97-158:`) updates pi according to this maximum-likelihood estimation formula:

```
pi[i] = (sum_j(nij[i,j]) + pc + root_state[i]) / (TINY + mu * sum_j(W[i,j]*Ti[j]) + sum(root_state) + pc)
```

The `root_state` vector appears in both numerator and denominator. It biases pi toward the composition observed at the tree root, providing a reasonable prior when substitution counts are sparse.

## How root_state is computed

Both v0 and v1 derive `root_state` from the root's marginal profile. The marginal profile is a 2D array of shape (L, n) where L is the alignment length and n is the number of states (4 for nucleotides). Each row is a probability distribution: `profile[pos, k]` is `P(state=k at position pos | all data)`.

These profiles come from Felsenstein's pruning algorithm (Felsenstein, 1981): a postorder traversal accumulates subtree likelihoods, then a preorder traversal combines them with outgroup information to produce the full marginal posterior.

To get `root_state`, the code takes the argmax of each row (most probable state) and counts how many positions map to each state.

## The problem: gap-only columns

An alignment column where all leaf sequences have gaps or ambiguous characters (N) carries no phylogenetic signal. No leaf constrains the ancestral state. The marginal profile at the root becomes uniform: `[1/n, 1/n, ..., 1/n]`.

For n=4 nucleotides, a gap-only column produces `profile[pos] = [0.25, 0.25, 0.25, 0.25]`.

The `argmax` of a uniform distribution is undefined mathematically. NumPy returns the lowest index (state 0). ndarray may return a different index due to floating-point differences in profile computation. Both are deterministic but implementation-dependent.

When v0 counts gap-only columns, it assigns them all to state 0 (or whichever state argmax selects). This adds an arbitrary bias to `root_state[0]` proportional to the number of gap-only columns in the alignment.

## v0: counts all positions

v0 counts every position unconditionally at `packages/legacy/treetime/treetime/treeanc.py:1608-1613:`:

```python
root_state = np.array([
    np.sum((self.tree.root.cseq == nuc) * self.data.multiplicity(mask=self.tree.root.mask))
    for nuc in self.gtr.alphabet
])
```

The `cseq` property returns the argmax sequence computed by `prof2seq()` at `packages/legacy/treetime/treetime/seq_utils.py:271:`:

```python
idx = tmp_profile.argmax(axis=1)
seq = gtr.alphabet[idx]
```

NumPy's `argmax` returns the first index among ties, so gap-only columns contribute to `root_state[0]` (state A for nucleotides).

## v1: skips uninformative positions

v1 checks whether each profile row has a dominant state before counting, at `packages/treetime/src/gtr/infer_gtr/dense.rs:144-152:`:

```rust
let uniform_threshold = 1.0 / n_states as f64 + 1e-10;
for row in root_profile.rows() {
    let max_val = row.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if max_val > uniform_threshold {
        let state = argmax_first(&row).ok_or_else(|| ...)?;
        counts[state] += 1.0;
    }
}
```

For n=4, the threshold is `0.25 + 1e-10`. A gap-only column with `max = 0.25` fails the test and is excluded. An informative column (e.g., `[0.95, 0.02, 0.02, 0.01]`) has `max = 0.95 > 0.25 + 1e-10` and is included.

The `1e-10` epsilon accounts for floating-point roundoff. Profiles are computed through matrix exponentiations and normalizations that may not produce exactly `1/n`.

## Scientific rationale

Gap-only columns contain no evolutionary information. Their presence depends on alignment quality and sequencing coverage, not on the substitution process. Including them in `root_state` adds a prior that:

1. Varies unpredictably across datasets based on alignment gaps
2. Is assigned to an arbitrary state based on argmax tie-breaking
3. Has no biological interpretation

v1 excludes these positions because they should not influence equilibrium frequency estimation. The remaining informative positions provide a biologically meaningful prior.

Both approaches produce valid GTR models. The difference is whether alignment artifacts (gap columns) influence the prior on pi.

## Golden master tests

The capture script at `packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture:` replicates v0's nij/Ti accumulation from `packages/legacy/treetime/treetime/treeanc.py:1556-1572:` but applies the same uninformative-position filter to root_state. It calls `GTR.infer()` at `packages/legacy/treetime/treetime/gtr.py:492-599:` directly rather than `tt.infer_gtr()`, which would include v0's unfiltered root_state computation.

Test cases in `packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs:` compare v1 output against the modified capture. All seven real datasets pass at 1e-6 tolerance. The tolerance was widened from 1e-7 to accommodate BLAS drift on mpox_clade_ii_20 (~200k positions).

## References

- Tavaré, S. (1986). Some probabilistic and statistical problems in the analysis of DNA sequences. _Lectures on Mathematics in the Life Sciences_, 17, 57-86.
- Felsenstein, J. (1981). Evolutionary trees from DNA sequences: A maximum likelihood approach. _Journal of Molecular Evolution_, 17(6), 368-376.
