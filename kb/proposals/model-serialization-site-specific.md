# Serialization format for site-specific model parameters

## Motivation

Candidate option R5 in the [auto-partitioning report](../reports/auto-partitioning.md) considers exporting inferred site-specific model parameters ($\pi^a$, $\mu^a$, $W$) in a reusable format, following IQ-TREE's two-stage PMSF workflow. The parameters require iterative approximate maximum-likelihood estimation over the tree but are cheaper to apply once fixed. Serializing them could enable:

- Reuse across analyses: estimate once, apply to rerooting, molecular clock, bootstrapping
- Cross-session persistence: avoid re-estimation when rerunning with different settings
- Inspection and debugging: verify inferred preferences and rates
- Interoperability: exchange profiles with IQ-TREE via `.sitefreq` export

TreeTime v0 has no serialization for site-specific GTR. The scalar `GTR` class has a text-based `__str__`/`from_file` round-trip, but `GTR_site_specific` has no save or load methods.

> **Status:** This proposal is not ready for implementation. Site-specific models have no production construction path, and the schema's scale normalization, invariant-safe loading, alphabet contract, and interoperability guarantees remain unapproved.

## Prior art

### IQ-TREE `.sitefreq`

Whitespace-delimited text, one line per alignment site. Each line contains a 1-based site index followed by exactly the model's canonical-state count of frequency values. IQ-TREE's writer left-aligns the index in a width-six field and uses width 15 for values, without a fixed 12-decimal contract. For amino acids, the column order is A R N D C Q E G H I L K M F P S T W Y V. The file stores $\pi$ only, not $\mu$ or $W$.

Example (amino acid, 20 states):

```
     1 0.050000000000 0.030000000000 0.040000000000 0.060000000000 0.010000000000 0.040000000000 0.070000000000 0.080000000000 0.020000000000 0.050000000000 0.090000000000 0.060000000000 0.020000000000 0.030000000000 0.050000000000 0.070000000000 0.060000000000 0.010000000000 0.030000000000 0.060000000000
     2 0.080000000000 0.020000000000 0.050000000000 0.040000000000 0.020000000000 0.030000000000 0.050000000000 0.100000000000 0.030000000000 0.060000000000 0.080000000000 0.050000000000 0.030000000000 0.040000000000 0.040000000000 0.060000000000 0.050000000000 0.010000000000 0.040000000000 0.070000000000
```

For TreeTime's canonical nucleotide model, each line would have four values in A C G T order; gap is an undetermined symbol outside the GTR state vector:

```
1      0.250000000000 0.200000000000 0.300000000000 0.250000000000
2      0.100000000000 0.100000000000 0.700000000000 0.100000000000
```

Reader validation: frequencies must be strictly positive and less than 1.0. Normalized if sum deviates from 1.0 (warning if > 1e-3). Low frequencies regularized via `convfreq()`. Sites with identical patterns and frequencies share a single profile entry.

Format details from [IQ-TREE source `readSiteStateFreq()`](https://github.com/iqtree/iqtree2/blob/main/alignment/alignment.cpp) (line 5620) and [`printSiteStateFreq()`](https://github.com/iqtree/iqtree2/blob/main/main/treetesting.cpp) (line 409).

### PhyloBayes `.siteprofiles`

`readpb_mpi -ss` exports posterior mean site-specific equilibrium frequencies from MCMC output to `<chainname>.siteprofiles`. Per-site amino acid frequency vectors averaged across MCMC samples. Only available under infinite mixture models (CAT, CAT-GTR). Requires `-s` flag during chain run to save model configurations. Related: `-rr` exports mean posterior relative exchangeabilities.

### RAxML-NG `.raxml.bestModel`

Per-partition model strings with inline optimized parameter values. One line per partition. Can be reused directly as `--model` input for subsequent runs.

Example (3 partitions):

```
GTR{0.200/1.000/2.000/4.000/7.000/1.000}+FU{0.280/0.220/0.230/0.270}+G4m{0.500}, NADH4 = 1-504
HKY{2.500}+FU{0.300/0.200/0.250/0.250}+G4m{1.200}, tRNA = 505-656
GTR{0.150/0.800/1.500/3.200/5.600/1.000}+FU{0.310/0.190/0.240/0.260}+I{0.150}+G4m{0.870}, NADH5 = 657-898
```

Format: `MODEL{rates}+FU{freqs}+G4m{alpha}` for user-supplied frequencies. `+FO` instead requests frequency optimization and takes no values. Gamma shape $\alpha$ is an among-site rate-heterogeneity parameter, not TreeTime's overall $\mu$. RAxML-NG provides no per-site frequency export in this format.

### BEAST2 XML state

Full model state serialized as XML. Parameters embedded in the model graph structure. Verbose but self-describing.

### TreeTime v0 scalar GTR text

`GTR.__str__()` produces human-readable multi-line text. `GTR.from_file()` parses it back via line-by-line string matching. Scalar GTR only, no site-specific variant exists.

Example (from v0 source code `gtr.py:145-170`):

```
Substitution rate (mu): 0.003421

Equilibrium frequencies (pi_i):
  A: 0.2983
  C: 0.1836
  G: 0.1963
  T: 0.3218

Symmetrized rates from j->i (W_ij):
     A  C  G  T
  A  0  1.2  2.5  0.8
  C  1.2  0  0.9  3.1
  G  2.5  0.9  0  0.7
  T  0.8  3.1  0.7  0
```

## Proposed format

### Primary: JSON

JSON stores the full site-specific model with metadata. TreeTime v1 already has JSON helpers (`json_write_str`, `json_write_file`, `json_read_str`) and ndarray serde support (`array1_as_vec`, `array2_as_vec`).

```json
{
  "format": "treetime-site-specific-gtr",
  "version": 1,
  "alphabet": ["A", "C", "G", "T"],
  "seq_len": 10000,
  "W": [[0, 1.2, 2.5, 0.8], ...],
  "pi": [[0.25, 0.30, ...], ...],
  "mu": [1.0, 0.5, 2.3, ...],
  "metadata": {
    "source_alignment": "data/flu/h3n2/200/aln.fasta",
    "source_tree": "data/flu/h3n2/200/tree.nwk",
    "inference_iterations": 15,
    "pseudocount": 1.0,
    "average_rate": 1.0
  }
}
```

Fields:

- `format`: format identifier for validation
- `version`: schema version for forward compatibility
- `alphabet`: state labels in column order
- `seq_len`: alignment length (redundant with array dimensions, for validation)
- `W`: shared exchangeability matrix (`n_states x n_states`, symmetric, zero diagonal)
- `pi`: per-site equilibrium frequencies (`n_states x seq_len`, columns sum to 1.0)
- `mu`: per-site rates (`seq_len`); the schema must define their weighted normalization and their scale relationship with $W$
- `metadata`: provenance and inference settings

### Secondary: IQ-TREE-compatible `.sitefreq` export

For interoperability, export $\pi$ in IQ-TREE's format:

```
     1 0.250000000000 0.300000000000 0.200000000000 0.250000000000
     2 0.100000000000 0.100000000000 0.700000000000 0.100000000000
```

This loses $\mu$, $W$, and metadata, but enables IQ-TREE to consume TreeTime's inferred profiles for ML tree search, and TreeTime to consume IQ-TREE's PMSF output.

### Scalar GTR

The same JSON format without per-site dimensions:

```json
{
  "format": "treetime-gtr",
  "version": 1,
  "alphabet": ["A", "C", "G", "T"],
  "W": [[0, 1.2, 2.5, 0.8], ...],
  "pi": [0.25, 0.30, 0.20, 0.25],
  "mu": 1.0
}
```

Replaces v0's text-based `__str__`/`from_file` with a structured format that is parseable without fragile line-by-line string matching.

## Design decisions

- JSON over binary (NPZ, HDF5): human-readable, debuggable, no additional dependencies. Per-site arrays for a 30K-site nucleotide alignment produce ~2 MB JSON (5 states x 30K sites x ~15 chars per float). Acceptable for TreeTime's typical dataset sizes
- JSON over YAML/TOML: TreeTime already has JSON infrastructure. Arrays of arrays map naturally to JSON. YAML/TOML handle nested numeric arrays poorly
- Separate `.sitefreq` export over embedding IQ-TREE format as primary: TreeTime needs $\mu$ and $W$ for its own reuse workflow. IQ-TREE's format is too limited as a primary format but valuable as an interoperability bridge
- `format` and `version` fields: enable schema evolution and validation without relying on file extension

## Impact

- Enables R5 (model serialization) from the auto-partitioning report
- Enables the two-stage PMSF-like workflow (section 7.5 of the report)
- Provides interoperability with IQ-TREE for site frequency profiles
- Replaces v0's fragile text parsing with structured JSON

## Validation plan

- Round-trip test: write a normalized model, reconstruct it through invariant-safe production constructors, and compare transition probabilities as well as fields
- IQ-TREE interoperability test: validate against a pinned IQ-TREE version and compare parsed profiles, rather than relying on manual acceptance
- Schema validation: reject files with wrong `format`, wrong `version`, mismatched dimensions
- Large dataset test: serialize and deserialize a 30K-site model, verify performance is acceptable

## Related

- [Auto-partitioning report](../reports/auto-partitioning.md) - design recommendation R5
- [Multi-partition config file proposal](config-file-multi-partition.md) - complementary I/O proposal
