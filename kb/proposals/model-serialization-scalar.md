# Serialization format for scalar GTR model parameters

## Motivation

TreeTime v1 infers GTR substitution model parameters ($W$, $\pi$, $\mu$) from sequence data but has no structured format to save them. Users who want to reuse an inferred model across analyses, inspect estimated parameters, or transfer models between tools must re-run inference each time. TreeTime v0 has a text-based `GTR.__str__()`/`GTR.from_file()` round-trip, but the format is fragile (line-by-line string matching, no metadata, no version info) and has no cross-tool compatibility.

A structured serialization format for scalar GTR enables:

- Reuse across analyses: infer model once, apply to subsequent runs via `--model-file`
- Cross-tool compatibility: import models from RAxML-NG or PAML, export for use in other tools
- Inspection and debugging: verify inferred parameters in a human-readable format
- Reproducibility: record exact model parameters alongside analysis results

## Prior art

Surveyed in the [substitution model serialization formats report](../reports/substitution-model-serialization-formats.md). Key formats for scalar GTR:

### PAML triangular format

De facto cross-tool standard for amino acid exchangeability matrices. Lower-triangular $W$ (190 values for 20 amino acids, or 6 values for 4 nucleotides) followed by equilibrium frequencies. Accepted by IQ-TREE, PhyML, RAxML-NG, and MEGA. Does not store $\mu$ (overall rate).

Example (nucleotide GTR, 4 states):

```
1.200000
2.500000 0.900000
0.800000 3.100000 0.700000

0.298300 0.183600 0.196300 0.321800
```

6 exchangeability values (AC, AG, AT, CG, CT, GT as lower triangle), blank line, 4 frequencies (A, C, G, T).

### RAxML-NG model string

Single-line format with inline parameter values. Stores $W$, $\pi$, $\mu$ (via rate heterogeneity parameters). Round-trip: RAxML-NG writes `.raxml.bestModel`, reads it back via `--model`.

Example:

```
GTR{0.200/1.000/2.000/4.000/7.000/1.000}+FO{0.280/0.220/0.230/0.270}+G4m{0.500}
```

### TreeTime v0 text

Multi-line human-readable format with labeled sections for $\mu$, $\pi$, $W$. Round-trip via `GTR.__str__()`/`GTR.from_file()`. Fragile parser depends on exact label text.

Example:

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

Stores the full scalar GTR model with metadata. Consistent with the [site-specific serialization proposal](model-serialization-site-specific.md), which uses the same structure with per-site array dimensions.

```json
{
  "format": "treetime-gtr",
  "version": 1,
  "alphabet": ["A", "C", "G", "T", "-"],
  "W": [
    [0.0, 1.2, 2.5, 0.8, 0.1],
    [1.2, 0.0, 0.9, 3.1, 0.1],
    [2.5, 0.9, 0.0, 0.7, 0.1],
    [0.8, 3.1, 0.7, 0.0, 0.1],
    [0.1, 0.1, 0.1, 0.1, 0.0]
  ],
  "pi": [0.2983, 0.1836, 0.1963, 0.3218, 0.0],
  "mu": 0.003421,
  "metadata": {
    "model_name": "GTR",
    "source_alignment": "data/flu/h3n2/200/aln.fasta",
    "inference_method": "marginal",
    "pseudocount": 1.0
  }
}
```

Fields:

- `format`: `"treetime-gtr"` (scalar) or `"treetime-site-specific-gtr"` (site-specific). Enables schema dispatch
- `version`: schema version for forward compatibility
- `alphabet`: state labels in array index order
- `W`: symmetric exchangeability matrix (`n_states x n_states`, zero diagonal)
- `pi`: equilibrium frequencies (`n_states`, sums to 1.0)
- `mu`: overall substitution rate (scalar)
- `metadata`: provenance and inference settings (optional fields, tool-specific)

### Secondary: PAML-compatible export

For cross-tool use, export $W$ and $\pi$ in PAML triangular format. This enables IQ-TREE, PhyML, and RAxML-NG to consume TreeTime's inferred models as custom matrices.

For nucleotides (4 states, 6 exchangeability values + 4 frequencies):

```
1.200000
2.500000 0.900000
0.800000 3.100000 0.700000

0.298300 0.183600 0.196300 0.321800
```

For nucleotides with gap state (5 states, 10 exchangeability values + 5 frequencies):

```
1.200000
2.500000 0.900000
0.800000 3.100000 0.700000
0.100000 0.100000 0.100000 0.100000

0.298300 0.183600 0.196300 0.321800 0.000000
```

PAML format does not store $\mu$. The rate scalar is absorbed into branch lengths in tools that consume this format.

### Secondary: RAxML-NG-compatible model string

For interoperability with RAxML-NG, export as a single-line model string:

```
GTR{1.200/2.500/0.800/0.900/3.100/0.700}+FO{0.298/0.184/0.196/0.322}
```

Rate order: AC/AG/AT/CG/CT/GT (upper triangle, row-major). This can be passed directly as `--model` input to RAxML-NG.

## Design decisions

- JSON as primary: consistent with site-specific proposal, human-readable, round-trip safe, supports metadata. TreeTime v1 already has JSON helpers (`json_write_str`, `json_read_str`, `JsonPretty`)
- `format` field discriminates scalar vs site-specific: a reader can dispatch on `"treetime-gtr"` vs `"treetime-site-specific-gtr"` without inspecting array dimensions
- PAML export for cross-tool use: the triangular format is the only format accepted by multiple tools. Worth supporting even though it loses $\mu$
- RAxML-NG string export for partition workflows: enables embedding TreeTime's inferred models in RAxML-NG partition files

## Impact

- Enables `--model-file` CLI flag for loading pre-inferred models
- Enables `--export-model` for saving inferred models after ancestral reconstruction or timetree
- Replaces v0's fragile text parser with structured JSON
- Provides PAML and RAxML-NG interoperability for scalar models
- Pairs with the [site-specific serialization proposal](model-serialization-site-specific.md) to cover both model variants

## Validation plan

- Round-trip test: write GTR to JSON, read back, compare all fields at 1e-15 tolerance
- PAML export test: write PAML triangular, verify manually in IQ-TREE with `-m <file>`
- RAxML-NG string test: write model string, verify RAxML-NG accepts it via `--model`
- Schema validation: reject files with wrong `format`, wrong `version`, negative frequencies, non-symmetric $W$
- v0 compatibility test: read v0's `GTR.from_file()` text format as a legacy import path

## Related

- [Serialization format for site-specific model parameters](model-serialization-site-specific.md) - sibling proposal for per-site models
- [Substitution model serialization formats report](../reports/substitution-model-serialization-formats.md) - format survey
- [Auto-partitioning report](../reports/auto-partitioning.md) - design recommendation R5
- [Multi-partition config file proposal](config-file-multi-partition.md) - per-partition model assignment references serialized models
