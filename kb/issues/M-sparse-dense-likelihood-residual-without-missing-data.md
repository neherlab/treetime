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
