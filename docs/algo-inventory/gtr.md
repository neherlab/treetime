# GTR Substitution Models

[Back to index](index.md)

## Substitution Models

### JC69 (Jukes-Cantor 1969)

| Property    | Value                                                                               |
| ----------- | ----------------------------------------------------------------------------------- |
| Type        | Well-known                                                                          |
| v1 Location | `packages/treetime/src/gtr/get_gtr.rs:84-102:`                                      |
| Reference   | Jukes & Cantor (1969). "Evolution of Protein Molecules." Academic Press, pp. 21-132 |

Simplest model: equal rates, equal frequencies.

### K80 (Kimura 2-Parameter 1980)

| Property    | Value                                           |
| ----------- | ----------------------------------------------- |
| Type        | Well-known                                      |
| v1 Location | `packages/treetime/src/gtr/get_gtr.rs:120-142:` |
| Reference   | Kimura (1980). J Mol Evol, 16(2):111-120        |
| Paper URL   | https://doi.org/10.1007/BF01731581              |

Distinguishes transitions from transversions.

### HKY85 (Hasegawa-Kishino-Yano 1985)

| Property    | Value                                             |
| ----------- | ------------------------------------------------- |
| Type        | Well-known                                        |
| v1 Location | `packages/treetime/src/gtr/get_gtr.rs:196-221:`   |
| Reference   | Hasegawa et al. (1985). J Mol Evol, 22(2):160-174 |
| Paper URL   | https://doi.org/10.1007/BF02101694                |

Combines K80 transition/transversion distinction with F81 unequal frequencies.

**Defect**: v1 lacks pi parameter (always uses uniform frequencies).

### TN93 (Tamura-Nei 1993)

| Property    | Value                                                 |
| ----------- | ----------------------------------------------------- |
| Type        | Well-known                                            |
| v1 Location | `packages/treetime/src/gtr/get_gtr.rs:386-422:`       |
| Reference   | Tamura & Nei (1993). Mol Biol Evol, 10(3):512-526     |
| Paper URL   | https://doi.org/10.1093/oxfordjournals.molbev.a040023 |

Two transition rate parameters (purine vs pyrimidine transitions).

### JTT92 (Jones-Taylor-Thornton 1992)

| Property    | Value                                           |
| ----------- | ----------------------------------------------- |
| Type        | Well-known (empirical AA model)                 |
| v1 Location | `packages/treetime/src/gtr/get_gtr.rs:280-350:` |
| Reference   | Jones et al. (1992). CABIOS, 8(3):275-282       |

20x20 empirical amino acid substitution matrix.

### Additional Models

- **F81** (`get_gtr.rs:156-178`): Unequal frequencies, equal rates
- **T92** (`get_gtr.rs:243-267`): GC-content parameterization

---

## Matrix Exponentiation

| Property    | Value                                                                                                       |
| ----------- | ----------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                  |
| v1 Location | `packages/treetime/src/gtr/gtr.rs:18-39:` (eig), `207-215:` (expQt)                                         |
| Reference   | Felsenstein (2004). "Inferring Phylogenies." Chapter 13                                                     |
| Reference   | Moler & Van Loan (2003). "Nineteen Dubious Ways to Compute the Matrix Exponential." SIAM Review, 45(1):3-49 |

Symmetrization trick: `S = sqrt(Pi) * Q * sqrt(Pi)^{-1}`, then `exp(Qt) = V * diag(exp(lambda_i * t)) * V^{-1}`.

---

## GTR Inference

| Property    | Value                                                |
| ----------- | ---------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                           |
| v1 Location | `packages/treetime/src/gtr/infer_gtr/mod.rs:70-127:` |
| v0 Location | `packages/legacy/treetime/treetime/gtr.py:491-599:`  |

Iterative coordinate descent: update W, normalize, update pi, update mu.

---

## Unimplemented

See [unimplemented.md](unimplemented.md) for full details:

- Site-specific GTR (`gtr_site_specific.py`)
- Random GTR generation
- File-based GTR loading

---

## File Index

| File                                       | Algorithms                              |
| ------------------------------------------ | --------------------------------------- |
| `packages/treetime/src/gtr/gtr.rs`         | GTR core, eigendecomposition, expQt     |
| `packages/treetime/src/gtr/get_gtr.rs`     | JC69, K80, F81, HKY85, T92, TN93, JTT92 |
| `packages/treetime/src/gtr/infer_gtr/*.rs` | GTR inference from mutation counts      |
