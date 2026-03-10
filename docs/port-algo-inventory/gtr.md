# GTR Substitution Models

[Back to index](../index.md)

## Substitution Models

### JC69 (Jukes-Cantor 1969)

| Property    | Value                                                                                                                             |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                        |
| v1 Location | `jc69()` (`#jc69`) in [`packages/treetime/src/gtr/get_gtr.rs#L161-L173`](../../../packages/treetime/src/gtr/get_gtr.rs#L161-L173) |
| Reference   | Jukes & Cantor (1969). "Evolution of Protein Molecules." Academic Press, pp. 21-132                                               |

Simplest model: equal rates, equal frequencies.

### K80 (Kimura 2-Parameter 1980)

| Property    | Value                                                                                                                           |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                      |
| v1 Location | `k80()` (`#k80`) in [`packages/treetime/src/gtr/get_gtr.rs#L200-L213`](../../../packages/treetime/src/gtr/get_gtr.rs#L200-L213) |
| Reference   | Kimura (1980). J Mol Evol, 16(2):111-120                                                                                        |
| Paper URL   | https://doi.org/10.1007/BF01731581                                                                                              |

Distinguishes transitions from transversions.

### HKY85 (Hasegawa-Kishino-Yano 1985)

| Property    | Value                                                                                                                               |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                          |
| v1 Location | `hky85()` (`#hky85`) in [`packages/treetime/src/gtr/get_gtr.rs#L280-L294`](../../../packages/treetime/src/gtr/get_gtr.rs#L280-L294) |
| Reference   | Hasegawa et al. (1985). J Mol Evol, 22(2):160-174                                                                                   |
| Paper URL   | https://doi.org/10.1007/BF02101694                                                                                                  |

Combines K80 transition/transversion distinction with F81 unequal frequencies. Accepts optional `pi` parameter for non-uniform equilibrium frequencies.

### TN93 (Tamura-Nei 1993)

| Property    | Value                                                                                                                             |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                        |
| v1 Location | `tn93()` (`#tn93`) in [`packages/treetime/src/gtr/get_gtr.rs#L466-L495`](../../../packages/treetime/src/gtr/get_gtr.rs#L466-L495) |
| Reference   | Tamura & Nei (1993). Mol Biol Evol, 10(3):512-526                                                                                 |
| Paper URL   | https://doi.org/10.1093/oxfordjournals.molbev.a040023                                                                             |

Two transition rate parameters (purine vs pyrimidine transitions).

### JTT92 (Jones-Taylor-Thornton 1992)

| Property    | Value                                                                                                                               |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (empirical AA model)                                                                                                     |
| v1 Location | `jtt92()` (`#jtt92`) in [`packages/treetime/src/gtr/get_gtr.rs#L360-L423`](../../../packages/treetime/src/gtr/get_gtr.rs#L360-L423) |
| Reference   | Jones et al. (1992). CABIOS, 8(3):275-282                                                                                           |

20x20 empirical amino acid substitution matrix.

### Additional Models

- **F81** `f81()` (`#f81`) in [`packages/treetime/src/gtr/get_gtr.rs#L237-L250`](../../../packages/treetime/src/gtr/get_gtr.rs#L237-L250): Unequal frequencies, equal rates. Accepts optional `pi` parameter.
- **T92** `t92()` (`#t92`) in [`packages/treetime/src/gtr/get_gtr.rs#L323-L340`](../../../packages/treetime/src/gtr/get_gtr.rs#L323-L340): GC-content parameterization

---

## Matrix Exponentiation

| Property    | Value                                                                                                                                                                                                                                                                      |
| ----------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                                                                                                                                                                 |
| v1 Location | `eig_single_site()` (`#eig_single_site`) in [`packages/treetime/src/gtr/gtr.rs#L72-L103`](../../../packages/treetime/src/gtr/gtr.rs#L72-L103), `expQt()` (`#expQt`) in [`packages/treetime/src/gtr/gtr.rs#L375-L385`](../../../packages/treetime/src/gtr/gtr.rs#L375-L385) |
| Reference   | Felsenstein (2004). "Inferring Phylogenies." Chapter 13                                                                                                                                                                                                                    |
| Reference   | Moler & Van Loan (2003). "Nineteen Dubious Ways to Compute the Matrix Exponential." SIAM Review, 45(1):3-49                                                                                                                                                                |

Symmetrization trick: `S = sqrt(Pi) * Q * sqrt(Pi)^{-1}`, then `exp(Qt) = V * diag(exp(lambda_i * t)) * V^{-1}`.

---

## GTR Inference

### Iterative Coordinate Descent

| Property    | Value                                                                                                                                                                 |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                                                                                                                                            |
| v1 Location | `infer_gtr_impl()` (`#infer_gtr_impl`) in [`packages/treetime/src/gtr/infer_gtr/common.rs#L97-L158`](../../../packages/treetime/src/gtr/infer_gtr/common.rs#L97-L158) |
| v0 Location | [`packages/legacy/treetime/treetime/gtr.py#L491-L599`](../../../packages/legacy/treetime/treetime/gtr.py#L491-L599)                                                   |

Iterative coordinate descent: update W, normalize, update pi, update mu. Shared by both sparse and dense inference paths via `MutationCounts` (`#MutationCounts`) struct.

### Sparse GTR Inference

| Property    | Value                                                                                                                                                                   |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                                                                                                                                              |
| v1 Location | `infer_gtr_sparse()` (`#infer_gtr_sparse`) in [`packages/treetime/src/gtr/infer_gtr/sparse.rs#L12-L21`](../../../packages/treetime/src/gtr/infer_gtr/sparse.rs#L12-L21) |

Counts mutations from Fitch reconstruction: integer substitution counts for `nij`, branch-length-weighted composition for `Ti`, root composition from consensus sequence.

### Dense GTR Inference

| Property    | Value                                                                                                                                                               |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                                                                                                                                          |
| v1 Location | `infer_gtr_dense()` (`#infer_gtr_dense`) in [`packages/treetime/src/gtr/infer_gtr/dense.rs#L21-L27`](../../../packages/treetime/src/gtr/infer_gtr/dense.rs#L21-L27) |

Uses fractional expected counts from branch joint distributions. Requires pre-populated profiles (two reconstruction passes when `--model=infer --dense=true`). Key helper functions:

- `get_branch_mutation_matrix()` (`#get_branch_mutation_matrix`) in [`packages/treetime/src/gtr/infer_gtr/dense.rs#L40-L71`](../../../packages/treetime/src/gtr/infer_gtr/dense.rs#L40-L71): Computes posterior `P(child=i, parent=j | site)` from edge messages and transition matrix.
- `accumulate_mutation_counts()` (`#accumulate_mutation_counts`) in [`packages/treetime/src/gtr/infer_gtr/dense.rs#L81-L119`](../../../packages/treetime/src/gtr/infer_gtr/dense.rs#L81-L119): Accumulates `nij` substitution counts and `Ti` time-in-state from branch joint distributions.
- `get_mutation_counts_dense()` (`#get_mutation_counts_dense`) in [`packages/treetime/src/gtr/infer_gtr/dense.rs#L127-L196`](../../../packages/treetime/src/gtr/infer_gtr/dense.rs#L127-L196): Iterates edges to build `MutationCounts` with `SUPERTINY_NUMBER` floor on `expQt` and branch length clamping.

---

## GTR Output

| Property    | Value                                                                                                                                             |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | I/O                                                                                                                                               |
| v1 Location | `write_gtr_json()` (`#write_gtr_json`) in [`packages/treetime/src/gtr/get_gtr.rs#L76-L80`](../../../packages/treetime/src/gtr/get_gtr.rs#L76-L80) |

Writes `gtr.json` with model type (named/inferred/custom), model name, mu, pi, and W. Parameters logged at info level via `log_gtr()` (`#log_gtr`).

---

## Unimplemented

See [unimplemented](../unimplemented/index.md) for full details:

- Site-specific GTR ([`packages/legacy/treetime/treetime/gtr_site_specific.py`](../../../packages/legacy/treetime/treetime/gtr_site_specific.py))
- Random GTR generation
- File-based GTR loading

---

## File Index

| File                                                                                                      | Algorithms                                               |
| --------------------------------------------------------------------------------------------------------- | -------------------------------------------------------- |
| [`packages/treetime/src/gtr/gtr.rs`](../../../packages/treetime/src/gtr/gtr.rs)                           | GTR core, eigendecomposition, `expQt()` (`#expQt`)       |
| [`packages/treetime/src/gtr/get_gtr.rs`](../../../packages/treetime/src/gtr/get_gtr.rs)                   | JC69, K80, F81, HKY85, T92, TN93, JTT92, GTR output JSON |
| [`packages/treetime/src/gtr/infer_gtr/common.rs`](../../../packages/treetime/src/gtr/infer_gtr/common.rs) | `MutationCounts`, `InferGtrOptions`, `infer_gtr_impl()`  |
| [`packages/treetime/src/gtr/infer_gtr/sparse.rs`](../../../packages/treetime/src/gtr/infer_gtr/sparse.rs) | Sparse GTR inference from Fitch mutation counts          |
| [`packages/treetime/src/gtr/infer_gtr/dense.rs`](../../../packages/treetime/src/gtr/infer_gtr/dense.rs)   | Dense GTR inference from branch joint distributions      |
