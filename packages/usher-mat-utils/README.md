# usher-mat-utils

Read and write UShER mutation-annotated tree (MAT) files in Protocol Buffer format.

## Overview

Provide utilities for working with phylogenetic trees produced by [UShER](https://github.com/yatisht/usher) (Ultrafast Sample placement on Existing tRees). UShER stores mutation-annotated trees in a custom protobuf format that encodes the tree structure (as Newick), per-node mutations, condensed nodes (collapsed identical sequences), and metadata.

## Supported Formats

Three protobuf schemas are supported:

- **parsimony** - Standard UShER MAT format with tree structure, mutations, condensed nodes, and clade annotations
- **mutation_detailed** - Extended format with detailed mutation encoding and sample placement data
- **taxodium** - Format used by [Taxodium](https://github.com/theosanderson/taxodium) for web-based tree visualization

## API

### Type Aliases

```rust
pub type UsherTree = parsimony::Data;
pub type UsherTreeNode = parsimony::CondensedNode;
pub type UsherMutation = parsimony::Mut;
pub type UsherMutationList = parsimony::MutationList;
pub type UsherMetadata = parsimony::NodeMetadata;
```

### Methods

```rust
impl UsherTree {
    /// Collect all mutation positions across the tree into a sorted set.
    fn get_all_positions(&self) -> BTreeSet<i32>;
}
```

### Functions

```rust
// Read from bytes or reader
fn usher_mat_pb_read_bytes(buf: impl Buf) -> Result<UsherTree, Report>;
fn usher_mat_pb_read(reader: impl Read) -> Result<UsherTree, Report>;

// Write to buffer or writer
fn usher_mat_pb_write_bytes(buf: &mut impl BufMut, tree: &UsherTree) -> Result<(), Report>;
fn usher_mat_pb_write(writer: &mut impl Write, tree: &UsherTree) -> Result<(), Report>;
```

## Usage

```rust
use usher_mat_utils::{usher_mat_pb_read, usher_mat_pb_write, UsherTree};
use std::fs::File;

// Read a MAT protobuf file
let file = File::open("tree.pb")?;
let tree: UsherTree = usher_mat_pb_read(file)?;

// Access tree data
println!("Newick: {}", tree.newick);
println!("Nodes with mutations: {}", tree.node_mutations.len());
println!("Condensed nodes: {}", tree.condensed_nodes.len());

// Write back to file
let mut output = File::create("output.pb")?;
usher_mat_pb_write(&mut output, &tree)?;
```

## Data Structure

The main `UsherTree` (`parsimony::Data`) contains:

- `newick` - Newick tree string (may include branch lengths)
- `node_mutations` - Per-node mutation lists in preorder traversal order
- `condensed_nodes` - Mapping of tree node names to collapsed identical sequences
- `metadata` - Per-node clade annotations in preorder traversal order

Mutations encode nucleotide changes with:

- `position` - Genomic position
- `ref_nuc`, `par_nuc`, `mut_nuc` - Reference, parent, and mutated nucleotides (encoded as 0:A, 1:C, 2:G, 3:T)

## Notes

All generated protobuf types derive `Serialize` and `Deserialize` for serde compatibility.

## References

- [UShER GitHub](https://github.com/yatisht/usher)
- [UShER Wiki](https://usher-wiki.readthedocs.io/)
