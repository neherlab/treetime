# treetime-cli

Command-line interface for Treetime - maximum-likelihood phylodynamic inference.

## Overview

Provide phylogenetic analysis tools including ancestral sequence reconstruction, molecular clock inference, and timetree estimation. This crate builds two binaries:

- `treetime` - main CLI for phylodynamic analysis
- `convert` - tree format conversion utility

## Installation

Build from the workspace root:

```bash
cargo build --release -p treetime-cli
```

## Commands

### treetime

```
treetime <COMMAND> [OPTIONS]
```

#### timetree

Estimate time trees from an initial tree topology, date constraints, and alignment (optional).

```bash
treetime timetree --tree tree.nwk --dates metadata.tsv --outdir output/ alignment.fasta
```

Key options:

- `--tree` - input tree file (Newick, Nexus, or Phylip)
- `--dates` - CSV/TSV file with node dates
- `--clock-rate` - fixed molecular clock rate
- `--coalescent` - coalescent time scale in years
- `--max-iter` - maximum optimization iterations (default: 2)
- `--outdir` - output directory

#### ancestral

Reconstruct ancestral sequences and map mutations to the tree.

```bash
treetime ancestral --tree tree.nwk --outdir output/ alignment.fasta
```

Key options:

- `--method-anc` - reconstruction method: `joint`, `marginal`, or `parsimony`
- `--model` - substitution model (e.g., `jc69`, `infer`)
- `--dense` - use dense probability representation

#### clock

Perform molecular clock analysis via root-to-tip regression.

```bash
treetime clock --tree tree.nwk --dates metadata.tsv --outdir output/
```

Key options:

- `--reroot` - rerooting strategy: `least-squares`, `min-dev`, `oldest`
- `--keep-root` - preserve current root position
- `--clock-filter` - outlier detection threshold (default: 3.0 IQD)
- `--covariation` - account for shared ancestry in regression

#### optimize

Optimize branch lengths and likelihood given aligned sequences.

```bash
treetime optimize --tree tree.nwk --aln alignment.fasta --outdir output/
```

Key options:

- `--model` - substitution model
- `--max-iter` - maximum iterations (default: 20)
- `--dp` - convergence threshold (default: 0.01)

#### prune

Remove branches from a phylogenetic tree.

```bash
treetime prune --tree tree.nwk --outdir output/ --prune-short 1e-6
```

Key options:

- `--prune-short` - remove branches below length threshold
- `--prune-empty` - remove branches without mutations (requires alignment)
- `--prune-nodes-list` - comma-separated node names to remove

#### homoplasy

Detect homoplasies (parallel mutations) in phylogenetic trees.

```bash
treetime homoplasy --tree tree.nwk --outdir output/ alignment.fasta
```

Key options:

- `--num-mut` - number of mutations to display (default: 10)
- `--detailed` - generate detailed report
- `--drms` - TSV file with drug resistance mutation info

#### mugration

Reconstruct discrete ancestral states (e.g., geographic location, host).

```bash
treetime mugration --tree tree.nwk --states traits.tsv --attribute country --outdir output/
```

Key options:

- `--attribute` - trait column to reconstruct
- `--states` - CSV/TSV file with discrete character states
- `--weights` - equilibrium state probabilities

#### completions

Generate shell completions.

```bash
treetime completions bash > ~/.local/share/bash-completion/treetime
```

Supported shells: bash, elvish, fish, fig, powershell, zsh

### convert

Convert between tree file formats.

```bash
convert input.nwk -o output.json -w auspice
```

Supported formats:

- `auspice` - Auspice JSON
- `newick` / `nexus` - standard tree formats
- `phyloxml` - PhyloXML
- `phyloxml-json` - PhyloXML as JSON
- `phylo-graph` - graph JSON
- `mat-pb` / `mat-json` - Usher mutation-annotated tree

Options:

- `-r` / `--input-format` - input format (auto-detected from extension)
- `-w` / `--output-format` - output format (auto-detected from extension)

## Common Options

All commands support:

- `-j` / `--jobs` - number of parallel threads
- `-v` / `--verbose` - increase verbosity (repeat for more: `-vv`, `-vvv`)
- `-q` / `--quiet` - decrease verbosity

## Input Formats

- Alignments: FASTA (plain or compressed: gz, bz2, xz, zstd)
- Trees: Newick, Nexus, Phylip
- Metadata: CSV or TSV with header row

## Workspace Dependencies

- `treetime` - core library (commands, models, algorithms)
- `treetime-graph` - graph data structure
- `treetime-io` - JSON/YAML serialization
- `treetime-utils` - error handling, compression, logging

## References

- Documentation: https://treetime.readthedocs.io/en/stable/
- Publication: https://academic.oup.com/ve/article/4/1/vex042/4794731
