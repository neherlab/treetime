# util-newick

Newick and Nexus phylogenetic tree parser and writer. Reads all major annotation dialects (BEAST, NHX, plain comments), eNewick networks, and Nexus container files. Writes in configurable styles.

## Data model

The central type is `NewickGraph` -- a flat adjacency list representing a tree (or DAG for eNewick networks).

```
NewickGraph
  nodes: Vec<NewickNodeData>    -- indexed by usize
  edges: Vec<NewickEdgeEntry>   -- indexed by usize
  root:  usize                  -- index into nodes
  rooted: Option<bool>          -- [&R] / [&U] marker
```

Nodes and edges live in flat `Vec`s. Everything references them by index. A node's `children` field holds edge indices (not node indices). To get a child node, go through the edge:

```rust
let node = &graph.nodes[some_node_idx];
for &edge_idx in &node.children {
    let edge = &graph.edges[edge_idx];
    let child_node = &graph.nodes[edge.child];
    // edge.data has branch_length, branch_attrs, raw_comments
    // child_node has name, node_attrs, raw_comments, hybrid, children
}
```

### What's on nodes vs edges

The split follows BEAST2's canonical grammar: annotations before `:` belong to the node, annotations after `:` belong to the edge.

```
    taxon[&node_stuff]:[&edge_stuff]0.05[&also_edge_stuff]
          ^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
          node_attrs           branch_attrs (merged)
```

**`NewickNodeData`** fields:

| Field          | Type                            | Content                                                                        |
| -------------- | ------------------------------- | ------------------------------------------------------------------------------ |
| `name`         | `Option<String>`                | taxon or clade label. `None` for anonymous nodes                               |
| `node_attrs`   | `BTreeMap<String, NewickValue>` | structured `[&k=v]` or `[&&NHX:k=v]` annotations before `:`                    |
| `raw_comments` | `Vec<String>`                   | plain `[text]` comments (no `&` prefix), preserved verbatim including brackets |
| `hybrid`       | `Option<NewickHybrid>`          | eNewick `#H1` / `##LGT2` marker, if present                                    |
| `children`     | `Vec<usize>`                    | edge indices into `graph.edges`, in parse order                                |

**`NewickEdgeData`** fields:

| Field           | Type                            | Content                                  |
| --------------- | ------------------------------- | ---------------------------------------- |
| `branch_length` | `Option<f64>`                   | the number after `:`. `None` when absent |
| `branch_attrs`  | `BTreeMap<String, NewickValue>` | structured annotations after `:`         |
| `raw_comments`  | `Vec<String>`                   | plain comments on the branch             |

### Annotation values

`NewickValue` has four variants. Which you get depends on the dialect:

**BEAST** `[&k=v]` -- typed parsing:

| Input                                   | Variant                   |
| --------------------------------------- | ------------------------- |
| `TRUE`, `FALSE` (case-insensitive)      | `Boolean(bool)`           |
| starts with digit or `-`, parses as f64 | `Number(f64)`             |
| `{1.0,2.0,3.0}`                         | `Array(Vec<NewickValue>)` |
| everything else                         | `String(String)`          |

**NHX** `[&&NHX:k=v]` -- all values are `String`. No type inference.

**Plain** `[text]` (no `&` prefix) -- not parsed into attrs, goes to `raw_comments` verbatim.

Dialect is detected per-comment automatically. A single tree can mix dialects on different nodes.

### Leaf vs internal

There's no enum distinguishing leaves from internal nodes. A leaf has `children.is_empty() == true`. An internal node has children.

### Finding nodes

The graph is a flat `Vec`, so finding nodes by name is a linear scan:

```rust
let node = graph.nodes.iter().find(|n| n.name.as_deref() == Some("taxon_A"));
```

For index-based lookup (when you need the index for edge traversal):

```rust
let idx = graph.nodes.iter().position(|n| n.name.as_deref() == Some("taxon_A")).unwrap();
```

### Parent lookup

Edges store both `parent` and `child` indices. To find a node's parent:

```rust
let parent_edge = graph.edges.iter().find(|e| e.child == node_idx);
if let Some(edge) = parent_edge {
    let parent_node = &graph.nodes[edge.parent];
}
```

The root node has no incoming edge.

## Parsing

```rust
use util_newick::{newick_from_string, newick_from_reader};

// From string
let graph = newick_from_string("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();

// From file
let file = std::fs::File::open("tree.nwk").unwrap();
let graph = newick_from_reader(std::io::BufReader::new(file)).unwrap();
```

Underscores in unquoted labels are preserved verbatim. The Felsenstein convention of converting underscores to spaces is not applied, as it breaks name matching against metadata.

The parser accepts:

- Standard Newick with optional names, branch lengths, quoted labels
- BEAST annotations `[&prob=0.95,hpd={1.0,2.0}]`
- NHX annotations `[&&NHX:S=human:B=90]`
- Plain comments `[any text]`
- Rooting markers `[&R]` / `[&U]` before the tree
- Empty branches `(,)` -- produces anonymous leaf nodes
- eNewick hybrid markers `#H1`, `##LGT2`

### Accessing parsed data

```rust
let graph = newick_from_string("(A[&prob=0.95]:0.1,B:0.2)root;").unwrap();

// Root
let root = &graph.nodes[graph.root];
assert_eq!(root.name.as_deref(), Some("root"));

// Iterate children of root
for &edge_idx in &root.children {
    let edge = &graph.edges[edge_idx];
    let child = &graph.nodes[edge.child];
    println!("child={} length={:?}", child.name.as_deref().unwrap_or("?"), edge.data.branch_length);
}

// Read an annotation
let a = graph.nodes.iter().find(|n| n.name.as_deref() == Some("A")).unwrap();
if let Some(util_newick::NewickValue::Number(prob)) = a.node_attrs.get("prob") {
    println!("posterior probability: {prob}");
}
```

## Writing

```rust
use util_newick::{newick_to_string, newick_to_writer, NewickWriteOptions, NwkStyle};

let opts = NewickWriteOptions::default(); // Beast style, full precision

// To string (returns Result for NHX reserved-char errors)
let nwk = newick_to_string(&graph, &opts).unwrap();

// To file
let mut file = std::fs::File::create("out.nwk").unwrap();
newick_to_writer(&mut file, &graph, &opts).unwrap();
```

### Output styles

Plain -- maximum compatibility, strips all annotations and comments:

```
(A:0.1,B:0.2)root;
```

Beast (default) -- BEAST2 canonical placement, node attrs before `:`, branch attrs between `:` and length:

```
(A[&prob=0.95]:[&rate=1.5]0.1,B:0.2)root;
```

Nhx -- NHX colon-separated format, all values as strings:

```
(A[&&NHX:prob=0.95]:0.1[&&NHX:rate=1.5],B:0.2)root;
```

Raw comments are emitted after structured blocks in Beast and Nhx styles. Plain drops them.

### Float precision

```rust
// Full precision (default)
NewickWriteOptions { significant_digits: None, decimal_digits: None, ..Default::default() }

// 3 significant digits
NewickWriteOptions { significant_digits: Some(3), ..Default::default() }
// -> 0.123

// 5 decimal places
NewickWriteOptions { decimal_digits: Some(5), significant_digits: None, ..Default::default() }
// -> 0.12346
```

## Building trees programmatically

```rust
use util_newick::*;
use std::collections::BTreeMap;

let mut graph = NewickGraph::new();

// Create nodes -- add_node returns the index
let root = graph.add_node(NewickNodeData::new().with_name("root"));
graph.root = root;

let leaf_a = graph.add_node(NewickNodeData::new().with_name("A"));
let leaf_b = graph.add_node(NewickNodeData::new().with_name("B"));

// Create edges -- add_edge returns the edge index and updates parent.children
graph.add_edge(root, leaf_a, NewickEdgeData::new().with_length(0.1));
graph.add_edge(root, leaf_b, NewickEdgeData::new().with_length(0.2));

// Add annotations
graph.nodes[leaf_a].node_attrs.insert("prob".to_owned(), NewickValue::Number(0.95));

// Rooting marker
graph.rooted = Some(true);

let nwk = newick_to_string(&graph, &NewickWriteOptions::default()).unwrap();
// -> [&R](A[&prob=0.95]:0.1,B:0.2)root;
```

### Building edges with annotations

```rust
let mut branch_attrs = BTreeMap::new();
branch_attrs.insert("rate".to_owned(), NewickValue::Number(1.5));

graph.add_edge(root, leaf_a, NewickEdgeData {
    branch_length: Some(0.1),
    branch_attrs,
    raw_comments: Vec::new(),
    is_acceptor: false,
});
```

## Nexus files

Nexus is a container around Newick strings. The crate handles parsing the TREES block, TRANSLATE tables, and multi-tree files.

```rust
use util_newick::{nexus_from_string, nexus_to_string, NewickWriteOptions};

let input = "#NEXUS\n\
    Begin Trees;\n\
      Translate\n\
        1 Homo_sapiens,\n\
        2 Pan_troglodytes\n\
      ;\n\
      Tree primates = [&R] (1:0.1,2:0.2);\n\
    End;\n";

let trees = nexus_from_string(input).unwrap();

// Each tree has a name and a NewickGraph
let tree = &trees[0];
assert_eq!(tree.name, "primates");
assert_eq!(tree.graph.rooted, Some(true));

// TRANSLATE tokens are resolved -- names are full taxon names, not integers
let names: Vec<_> = tree.graph.nodes.iter()
    .filter_map(|n| n.name.as_deref())
    .collect();
assert!(names.contains(&"Homo_sapiens"));

// Write back to Nexus
let output = nexus_to_string(&trees, &NewickWriteOptions::default()).unwrap();
```

The Nexus parser:

- Handles case-insensitive keywords (`begin`, `Begin`, `BEGIN`)
- Resolves TRANSLATE integer-to-name mappings before Newick parsing
- Ignores unknown blocks (FigTree display settings, Data blocks, etc.)
- Supports multiple `Tree` entries per TREES block

The Nexus writer emits a Taxa block (sorted, deduplicated leaf names) and a Trees block.

## eNewick networks

eNewick encodes phylogenetic networks (DAGs with hybrid/reticulate nodes) by duplicating hybrid nodes with `#` markers:

```
((A,(B)x#H1)c,(x#H1,C)d);
```

The parser detects `#TypeN` patterns in labels, creates a single graph node for each unique `(type, index)` pair, and connects multiple parent edges to it. The resulting graph is a DAG, not a tree.

```rust
let g = newick_from_string("((A,(B)x#H1)c,(x#H1,C)d);").unwrap();

// x#H1 appears as one node with two parent edges
let hybrid = g.nodes.iter().find(|n| n.hybrid.is_some()).unwrap();
let h = hybrid.hybrid.as_ref().unwrap();
assert_eq!(h.kind.as_deref(), Some("H"));  // H, LGT, R, or None
assert_eq!(h.index, 1);
```

The writer reproduces hybrid markers: first visit writes the full subtree, subsequent visits write `name#TypeN` with branch length only.

Reticulation types: `H` (hybridization), `LGT` (lateral gene transfer), `R` (recombination), or any custom string. `##` (double hash) marks the acceptor edge (main lineage).

## Equality

`PartialEq` on `NewickGraph` is **order-insensitive**: `(A,B)` equals `(B,A)`. Children are sorted by a canonical key `(name, child_count, branch_length)` before comparison. This is the right default for comparing trees parsed from different sources.

For **order-sensitive** comparison (when child ordering matters, e.g. verifying exact serialization round-trips):

```rust
assert!(g1.eq_ordered(&g2));
```

Both comparisons check: topology, node names, branch lengths (bitwise f64 equality), structured annotations, raw comments, hybrid markers, and rooting.

## Round-trip guarantees

- `parse(write(g, Beast))` == `g` for any valid graph (order-insensitive)
- `write(parse(write(g)))` == `write(g)` (string-level idempotence)
- Plain style is lossy by design: annotations and raw comments are dropped
- NHX converts typed values to strings: `Number(0.95)` becomes `String("0.95")`
- Nexus TRANSLATE is transparent: parsed names are always resolved, never integer tokens
