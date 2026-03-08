# GTR Substitution Models

[Back to index](../index.md)

## Substitution Models

### JC69 (Jukes-Cantor 1969)

| Property    | Value                                                                                                                           |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                      |
| v1 Location | `jc69()` (`#jc69`) in [`packages/treetime/src/gtr/get_gtr.rs#L84-L102`](../../../packages/treetime/src/gtr/get_gtr.rs#L84-L102) |
| Reference   | Jukes & Cantor (1969). "Evolution of Protein Molecules." Academic Press, pp. 21-132                                             |

Simplest model: equal rates, equal frequencies.

### K80 (Kimura 2-Parameter 1980)

| Property    | Value                                                                                                                           |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                      |
| v1 Location | `k80()` (`#k80`) in [`packages/treetime/src/gtr/get_gtr.rs#L120-L142`](../../../packages/treetime/src/gtr/get_gtr.rs#L120-L142) |
| Reference   | Kimura (1980). J Mol Evol, 16(2):111-120                                                                                        |
| Paper URL   | https://doi.org/10.1007/BF01731581                                                                                              |

Distinguishes transitions from transversions.

### HKY85 (Hasegawa-Kishino-Yano 1985)

| Property    | Value                                                                                                                               |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                          |
| v1 Location | `hky85()` (`#hky85`) in [`packages/treetime/src/gtr/get_gtr.rs#L196-L221`](../../../packages/treetime/src/gtr/get_gtr.rs#L196-L221) |
| Reference   | Hasegawa et al. (1985). J Mol Evol, 22(2):160-174                                                                                   |
| Paper URL   | https://doi.org/10.1007/BF02101694                                                                                                  |

Combines K80 transition/transversion distinction with F81 unequal frequencies.

**Defect**: v1 lacks pi parameter (always uses uniform frequencies).

### TN93 (Tamura-Nei 1993)

| Property    | Value                                                                                                                             |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known                                                                                                                        |
| v1 Location | `tn93()` (`#tn93`) in [`packages/treetime/src/gtr/get_gtr.rs#L386-L422`](../../../packages/treetime/src/gtr/get_gtr.rs#L386-L422) |
| Reference   | Tamura & Nei (1993). Mol Biol Evol, 10(3):512-526                                                                                 |
| Paper URL   | https://doi.org/10.1093/oxfordjournals.molbev.a040023                                                                             |

Two transition rate parameters (purine vs pyrimidine transitions).

### JTT92 (Jones-Taylor-Thornton 1992)

| Property    | Value                                                                                                                               |
| ----------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Well-known (empirical AA model)                                                                                                     |
| v1 Location | `jtt92()` (`#jtt92`) in [`packages/treetime/src/gtr/get_gtr.rs#L280-L350`](../../../packages/treetime/src/gtr/get_gtr.rs#L280-L350) |
| Reference   | Jones et al. (1992). CABIOS, 8(3):275-282                                                                                           |

20x20 empirical amino acid substitution matrix.

### Additional Models

- **F81** `f81()` (`#f81`) in [`packages/treetime/src/gtr/get_gtr.rs#L156-L178`](../../../packages/treetime/src/gtr/get_gtr.rs#L156-L178): Unequal frequencies, equal rates
- **T92** `t92()` (`#t92`) in [`packages/treetime/src/gtr/get_gtr.rs#L243-L267`](../../../packages/treetime/src/gtr/get_gtr.rs#L243-L267): GC-content parameterization

---

## Matrix Exponentiation

| Property    | Value                                                                                                                                                                                                                                                                    |
| ----------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Type        | Well-known                                                                                                                                                                                                                                                               |
| v1 Location | `eig_single_site()` (`#eig_single_site`) in [`packages/treetime/src/gtr/gtr.rs#L18-L39`](../../../packages/treetime/src/gtr/gtr.rs#L18-L39), `expQt()` (`#expQt`) in [`packages/treetime/src/gtr/gtr.rs#L207-L215`](../../../packages/treetime/src/gtr/gtr.rs#L207-L215) |
| Reference   | Felsenstein (2004). "Inferring Phylogenies." Chapter 13                                                                                                                                                                                                                  |
| Reference   | Moler & Van Loan (2003). "Nineteen Dubious Ways to Compute the Matrix Exponential." SIAM Review, 45(1):3-49                                                                                                                                                              |

Symmetrization trick: `S = sqrt(Pi) * Q * sqrt(Pi)^{-1}`, then `exp(Qt) = V * diag(exp(lambda_i * t)) * V^{-1}`.

---

## GTR Inference

| Property    | Value                                                                                                                                                           |
| ----------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| Type        | Custom (TreeTime-specific)                                                                                                                                      |
| v1 Location | `infer_gtr_impl()` (`#infer_gtr_impl`) in [`packages/treetime/src/gtr/infer_gtr/mod.rs#L70-L127`](../../../packages/treetime/src/gtr/infer_gtr/mod.rs#L70-L127) |
| v0 Location | [`packages/legacy/treetime/treetime/gtr.py#L491-L599`](../../../packages/legacy/treetime/treetime/gtr.py#L491-L599)                                             |

Iterative coordinate descent: update W, normalize, update pi, update mu.

---

## Unimplemented

See [unimplemented](../unimplemented/index.md) for full details:

- Site-specific GTR ([`packages/legacy/treetime/treetime/gtr_site_specific.py`](../../../packages/legacy/treetime/treetime/gtr_site_specific.py))
- Random GTR generation
- File-based GTR loading

---

## File Index

| File                                                                                                | Algorithms                                         |
| --------------------------------------------------------------------------------------------------- | -------------------------------------------------- |
| [`packages/treetime/src/gtr/gtr.rs`](../../../packages/treetime/src/gtr/gtr.rs)                     | GTR core, eigendecomposition, `expQt()` (`#expQt`) |
| [`packages/treetime/src/gtr/get_gtr.rs`](../../../packages/treetime/src/gtr/get_gtr.rs)             | JC69, K80, F81, HKY85, T92, TN93, JTT92            |
| [`packages/treetime/src/gtr/infer_gtr/mod.rs`](../../../packages/treetime/src/gtr/infer_gtr/mod.rs) | GTR inference from mutation counts                 |
