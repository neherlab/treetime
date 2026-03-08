# Uninformative root_state filtering

| Property      | Value                                                                                                                                                                                    |
| ------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type          | Intentional deviation from v0                                                                                                                                                            |
| v1 Location   | `get_mutation_counts_dense()` (`#get_mutation_counts_dense`) in [`packages/treetime/src/gtr/infer_gtr/dense.rs#L135-L149`](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L135-L149) |
| v0 Location   | [`packages/legacy/treetime/treetime/treeanc.py#L1608-L1613`](../../packages/legacy/treetime/treetime/treeanc.py#L1608-L1613)                                                             |
| Affects       | `root_state` vector, and downstream `pi` (equilibrium frequencies), `W` (exchangeability), `mu` (rate)                                                                                   |
| Datasets      | Any alignment with gap-only columns. Measured: lassa_L_50 (~29% W shift), tb_20 (~1e-5), mpox_clade_ii_20                                                                                |
| Golden master | Capture script updated to match v1. Fixtures regenerated. All 7 real datasets enabled.                                                                                                   |

##What root_state does

Dense GTR inference computes `root_state[k]`: the number of alignment positions
in state k at the tree root. This vector enters the iterative GTR solver
(`infer_gtr_impl()` (`#infer_gtr_impl`) in
[`packages/treetime/src/gtr/infer_gtr/common.rs#L120-L121`](../../packages/treetime/src/gtr/infer_gtr/common.rs#L120-L121),
v0 equivalent at
[`packages/legacy/treetime/treetime/gtr.py#L571-L574`](../../packages/legacy/treetime/treetime/gtr.py#L571-L574))
as a prior on equilibrium frequencies:

```
pi[i] = (sum_j nij[i,j] + pc + root_state[i]) / (TINY + mu * sum_j W[i,j]*Ti[j] + sum(root_state) + pc)
```

Both v0 and v1 compute `root_state` by taking `argmax` of the root's marginal
profile at each alignment position, then counting how many positions map to each
state.

##The deviation

Gap-only alignment columns (all leaves are gap or N) produce a uniform marginal
profile `[1/n, 1/n, ..., 1/n]` at the root. `argmax` on a uniform row returns
an arbitrary state index. The specific index depends on BLAS implementation:
NumPy (v0) and ndarray (v1) break ties differently due to floating-point
differences in matrix exponentiation.

For datasets with few gap-only columns, the effect is negligible. For lassa_L_50
(557 gap-only columns out of ~10k), hundreds of counts shift to one arbitrary
state, causing ~29% shift in `W`.

Whether including uninformative positions is correct is debatable. v0 includes
them, treating every position equally regardless of signal content. v1 excludes
them, on the basis that gap-only positions carry zero phylogenetic signal and
their contribution to `root_state` is determined entirely by tie-breaking
behavior rather than data.

Both approaches produce valid GTR models. The scientific argument for exclusion:
uninformative positions add an arbitrary prior to `pi` proportional to the number
of gap-only columns, which varies across datasets and depends on alignment
quality rather than evolutionary signal.

##v0: counts all positions

v0 counts all positions unconditionally via the consensus sequence `cseq`
([`packages/legacy/treetime/treetime/treeanc.py#L1608-L1613`](../../packages/legacy/treetime/treetime/treeanc.py#L1608-L1613)):

```python
root_state = np.array([
    np.sum((self.tree.root.cseq == nuc) * self.data.multiplicity(mask=self.tree.root.mask))
    for nuc in self.gtr.alphabet
])
```

`cseq` is derived from `argmax(marginal_profile, axis=1)` via `prof2seq()`
([`packages/legacy/treetime/treetime/seq_utils.py#L271`](../../packages/legacy/treetime/treetime/seq_utils.py#L271)),
so uniform rows produce `cseq[pos] = alphabet[0]` (NumPy argmax returns first
maximum).

##v1: skips uninformative positions

v1 skips positions where `max(profile_row) <= 1/n + 1e-10`
([`packages/treetime/src/gtr/infer_gtr/dense.rs#L139-L147`](../../packages/treetime/src/gtr/infer_gtr/dense.rs#L139-L147)):

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

For n=4 (nucleotides), threshold = 0.25 + 1e-10. A gap-only column with profile
`[0.25, 0.25, 0.25, 0.25]` has max = 0.25, not above threshold, excluded. An
informative column (e.g. `[0.95, 0.02, 0.02, 0.01]`) has max = 0.95, included.

##Golden master tests

The capture script
([`packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture`](../../packages/treetime/src/gtr/infer_gtr/__tests__/__fixtures__/gm_infer_gtr_dense_capture))
replicates v0's internal nij/Ti accumulation
([`packages/legacy/treetime/treetime/treeanc.py#L1556-L1572`](../../packages/legacy/treetime/treetime/treeanc.py#L1556-L1572))
with the same uninformative-position filter applied to root_state, then calls
`GTR.infer()`
([`packages/legacy/treetime/treetime/gtr.py#L492-L599`](../../packages/legacy/treetime/treetime/gtr.py#L492-L599))
directly. It does not call `tt.infer_gtr()` because v0's `infer_gtr()` includes
the unfiltered root_state computation.

Test cases
([`packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs`](../../packages/treetime/src/gtr/infer_gtr/__tests__/test_gm_infer_gtr_dense.rs))
re-enabled after the fix:

- `lassa_L_50`: was ~29% W shift, now passes at 1e-6
- `mpox_clade_ii_20`: now passes at 1e-6
- `tb_20`: now passes at 1e-6

Real-dataset tolerance widened from 1e-7 to 1e-6 to accommodate mpox_clade_ii_20
(~200k positions), where BLAS drift between NumPy and ndarray accumulates to
~2.3e-7.

##v0 status

Fixing or modifying v0 is out of scope for this project. The golden master
capture script applies the same filter so that oracle outputs and v1 outputs
agree on the filtered algorithm.
