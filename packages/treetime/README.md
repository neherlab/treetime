# treetime

Core library for phylogenetic analysis via ancestral sequence reconstruction and molecular clock inference. Scale branch lengths to calendar time by aligning tips with sampling dates and internal nodes with divergence times.

## Features

- Ancestral sequence reconstruction using marginal likelihood methods
- Molecular clock analysis with root-to-tip regression
- Timetree inference with iterative optimization
- Support for dense and sparse sequence representations
- General Time Reversible (GTR) substitution models
- Probability distributions for temporal constraints
- Coalescent-based population dynamics modeling

## Module Overview

- `alphabet` - Sequence alphabet definitions (nucleotides, amino acids)
- `cli` - Command-line interface utilities
- `commands` - Analysis commands:
  - `ancestral` - Ancestral sequence reconstruction
  - `clock` - Molecular clock regression and rerooting
  - `homoplasy` - Homoplasy detection
  - `mugration` - Discrete trait (migration) analysis
  - `optimize` - Node position optimization
  - `prune` - Branch pruning utilities
  - `timetree` - Full timetree inference pipeline
- `constants` - Global constants
- `distribution` - Probability distributions and operations (convolution, multiplication, scaling)
- `graph` - Graph tests (implementation in `treetime-graph`)
- `gtr` - GTR substitution model implementation and inference
- `io` - File format parsers and writers:
  - FASTA, Newick, NEXUS, PhyloXML
  - Auspice JSON, Graphviz DOT
  - CSV/TSV for dates and metadata
- `representation` - Data structures for partitions:
  - Parsimony and marginal likelihood partitions
  - Dense and sparse sequence storage
  - Timetree node/edge data
  - Primitives (`Seq`, `AsciiChar`, `BitSet128`, `StateSet`) in `treetime-primitives`
- `seq` - Sequence utilities (mutations, indels, composition)

## Dependencies

This crate integrates with workspace crates:

- `treetime-primitives` - Primitive types (`Seq`, `AsciiChar`, `BitSet128`, `StateSet`)
- `treetime-graph` - Generic graph structure (nodes, edges, traversal)
- `treetime-ops` - Numerical convolution and multiplication algorithms
- `treetime-grid` - Grid utilities for discretized distributions
- `treetime-analytical` - Analytical distribution formulas
- `treetime-io` - Low-level I/O utilities
- `treetime-utils` - Shared error handling and helpers

External dependencies include `ndarray` for numerical operations, `petgraph` for graph algorithms, and `argmin` for optimization.
