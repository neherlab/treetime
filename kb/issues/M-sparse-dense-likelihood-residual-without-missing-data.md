# Sparse and dense log-likelihoods differ on alignments without gaps or `N`

On the same tree with the same branch lengths, the sparse marginal log-likelihood is lower than the dense one on some datasets that keep a difference after every column with a gap or `N` is removed, so the difference does not come from missing data.

| Dataset | Alignment | Sparse minus dense |
|---|---|---|
| `rsv/a/20` | original | -0.62 |
| `rsv/a/20` | columns with gaps or `N` removed | -0.45 |
| `dengue/20` | original | -0.13 |
| `dengue/20` | columns with gaps or `N` removed | -0.12 |
| `lassa/L/20`, `ebola/20` | columns with gaps or `N` removed | 0 |
| `sc2/2844` | original, final `optimize` tree | -5.0 |

The differences move `optimize` tip lengths by up to 1.3% between the two representations on `dengue/20` and leave node dates of `rsv/a/20` and `lassa/L/20` slightly apart. `rsv/a/20` and `dengue/20` still contain IUPAC ambiguity codes after the column filter.

> [!IMPORTANT]
> **Investigation required.** Two causes fit the data: the folding approximation of [kb/decisions/ancestral-marginal-sparse-parsimony-chain.md](../decisions/ancestral-marginal-sparse-parsimony-chain.md) (posterior above 1 - 1e-4), and the handling of IUPAC codes in sparse messages. Compare sparse and dense per site on `rsv/a/20` after the column filter, once with the folding threshold disabled and once with the IUPAC codes replaced by `N`, to tell them apart.

## Related issues

- [M-sparse-non-root-posterior-folds-sites-where-a-child-disagrees.md](M-sparse-non-root-posterior-folds-sites-where-a-child-disagrees.md): below the root, the forward pass also folds substitution sites; on `data/rsv/a/20`, the exact rule of that issue moves the `optimize` branch lengths by at most 3.9e-7, while sparse and dense still differ by up to 2.0e-5, so that fold does not explain the residual there

## Proposed implementation plan

The table values come from earlier runs and were not measured again on the current code. The plan starts with the investigation of the block above, then fixes the cause it finds.

1. Write a test-only comparison in `packages/treetime/src/partition/marginal/__tests__/` (or reuse the sparse-dense harness of `packages/treetime/src/optimize/__tests__/test_topology_cleanup_sparse_dense.rs`) that computes per-site log-likelihoods of sparse and dense marginal inference on the same tree and branch lengths, and lists the sites whose difference exceeds `1e-9`
2. Run it on `data/rsv/a/20` after the column filter in two variants: folding switched off (the fold sites of `packages/treetime/src/partition/marginal/sparse/mutations.rs` and `packages/treetime/src/partition/marginal/shared/pass.rs`), and IUPAC codes replaced by `N`. The variant that removes the difference names the cause
3. Fold cause: the folding threshold is an approved approximation ([kb/decisions/ancestral-marginal-sparse-parsimony-chain.md](../decisions/ancestral-marginal-sparse-parsimony-chain.md)); report the per-site error to the user before any change of the decision. IUPAC cause: make sparse messages use the same ambiguity profile as dense at the sites with IUPAC codes
4. Regression test: dense versus sparse total log-likelihood on a small hand-made alignment with IUPAC codes and no gaps, equal within `1e-10`. Oracle: dense inference (independent path)

Validation: `just test-rs`, `just lint-rs`, then `just smoke` on `optimize` for `rsv/a/20` and `dengue/20` per `cli-verification.instructions.md`.
