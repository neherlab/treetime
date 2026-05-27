# Per-gene amino acid reconstruction for node data JSON

## Motivation

Node data JSON from `augur ancestral --genes` includes per-gene AA mutations (`aa_muts`), root AA sequences (`aa_sequences`), gene coordinate annotations, and AA reference sequences. This requires running ancestral reconstruction once per gene on pre-translated AA alignments. v1's reconstruction pipeline is alphabet-agnostic (works with AA alphabet and JTT92 model) but lacks the per-gene orchestration layer.

## What augur does

`augur/ancestral.py:580-685` `reconstruct_translations()`:

1. Load gene features from GenBank/GFF annotation file (CDS coordinates)
2. For each gene:
   a. Load pre-translated AA alignment FASTA
   b. Optionally load AA root sequence (from file, or translate from nuc reference)
   c. Call `run_ancestral(T, aa_aln, alphabet='aa')` - the same reconstruction function as nuc
   d. Merge per-gene mutations into the nuc JSON:
   - `nodes.<name>.aa_muts.<gene>` = gene's mutation list
   - `nodes.<root>.aa_sequences.<gene>` = root AA sequence (root only)
   - `reference.<gene>` = AA root/reference sequence
   - `annotations.<gene>` = CDS coordinates from GenBank/GFF feature

This is augur-side orchestration. TreeTime v0 does not do this - it just provides `TreeAnc` as a library. The v0 standalone CLI's `--aa` flag runs a single reconstruction with AA alphabet, no per-gene loop.

## What v1 has

- AA alphabet (`AlphabetName::Aa`, `AlphabetName::AaNoStop`)
- JTT92 substitution model
- Alphabet-agnostic reconstruction pipeline (marginal + parsimony)
- `--aa` CLI flag (parsed, not wired)

## What v1 lacks

- `--aa` flag not wired to select AA alphabet in the ancestral pipeline
- No GenBank/GFF gene annotation parsing
- No per-gene partition support (loop over genes, run reconstruction per gene)
- No `--translations` / `--annotation` / `--aa-root-sequence` CLI flags wired
- No merge of per-gene results into a combined output structure

## Implementation

1. Wire `--aa` flag to select AA alphabet and JTT92 (or inferred) model for single-partition reconstruction
2. Add GenBank/GFF feature parsing (BioPython equivalent for extracting CDS coordinates)
3. Add `--genes`, `--translations`, `--annotation`, `--aa-root-sequence` CLI flags
4. Implement per-gene loop: for each gene, create a partition with AA alphabet, run reconstruction, collect mutations
5. Merge AA results into the node data JSON output alongside nuc results
6. Add `aa_muts` and `aa_sequences` fields to the node data JSON serde types

No new reconstruction algorithm is needed. The core algorithm is identical to nucleotide - different alphabet size (20 vs 5) and model (JTT92 vs JC69/inferred).

## Validation

- Compare v1 `aa_muts` output against augur's `ancestral --genes` on test datasets
- Verify gene annotations match augur's `genome_features_to_auspice_annotation()` output
- Verify `aa_sequences` on root matches augur's output
