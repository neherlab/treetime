# treetime

Core library for phylogenetic analysis via ancestral sequence reconstruction and molecular clock inference. Scales branch lengths to calendar time by aligning tips with sampling dates and internal nodes with divergence times.

## Features

- Ancestral sequence reconstruction using marginal likelihood (dense and sparse) and parsimony methods
- Molecular clock analysis with root-to-tip regression and automatic rerooting
- Timetree inference with iterative optimization and coalescent modeling
- General Time Reversible (GTR) substitution model fitting
- Dense and sparse sequence representations
- Discrete trait (mugration) analysis

## Module Overview

- `alphabet` - Sequence alphabet definitions (nucleotides, amino acids) and configuration
- `cli` - CLI rendering utilities (root-to-tip charts)
- `commands` - Analysis pipelines:
  - `ancestral` - Ancestral sequence reconstruction (Fitch parsimony and marginal likelihood)
  - `clock` - Molecular clock regression, date constraints, rerooting
  - `homoplasy` - Homoplasy detection
  - `mugration` - Discrete trait (migration) analysis
  - `optimize` - Branch length optimization (dense and sparse)
  - `prune` - Branch pruning and collapsing
  - `timetree` - Full timetree inference: coalescent modeling, convergence tracking, forward/backward passes, clock filtering, polytomy resolution, relaxed clock
- `constants` - Numerical constants (branch length cutoffs)
- `graph` - Integration tests for graph operations (implementation in `treetime-graph`)
- `gtr` - GTR substitution model construction and inference
- `hacks` - Workarounds (branch length fixes)
- `io` - Integration tests for I/O (implementation in `treetime-io`)
- `representation` - Partition system and node/edge payloads:
  - `payload/` - Per-node and per-edge data (ancestral, timetree, dense, sparse)
  - `partition/` - Partition types (Fitch, marginal dense/sparse, timetree, likelihood) and traits
  - `algo/` - Inference algorithms (dense reconstruction)
- `seq` - Sequence utilities (mutations, indels, composition, divergence)

## Workspace Dependencies

- `treetime-primitives` - Primitive types (`Seq`, `AsciiChar`, `BitSet128`, `StateSet`)
- `treetime-graph` - Generic graph structure (nodes, edges, traversal)
- `treetime-distribution` - Probability distributions and operations (convolution, multiplication, scaling)
- `treetime-ops` - Numerical array operations (convolution, multiplication)
- `treetime-grid` - Grid interpolation for discretized distributions
- `treetime-analytical` - Analytical distribution formulas (gaussian, exponential)
- `treetime-io` - File format I/O (FASTA, Newick, NEXUS, PhyloXML, Auspice JSON, DOT, CSV/TSV)
- `treetime-validation` - Validation utilities
- `treetime-utils` - Error macros, compression, datetime helpers
- `phyloxml` - PhyloXML parsing
- `usher-mat-utils` - Usher MAT utilities

Key external dependencies: `ndarray` (numerical operations), `argmin` (optimization), `petgraph` (graph algorithms), `rayon` (parallelism).
