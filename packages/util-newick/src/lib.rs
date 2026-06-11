//! Newick and Nexus phylogenetic tree parser and writer.
//!
//! Reads and writes Newick trees with automatic dialect detection (BEAST `[&k=v]`,
//! NHX `[&&NHX:k=v]`, plain comments), eNewick network markers (`#H1`, `##LGT2`),
//! and Nexus container files (TRANSLATE tables, multi-tree blocks).
//!
//! # Parse a Newick string
//!
//! ```
//! use util_newick::newick_from_string;
//!
//! let graph = newick_from_string("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;").unwrap();
//!
//! // Traverse nodes
//! let root = &graph.nodes[graph.root];
//! assert_eq!(root.name.as_deref(), Some("F"));
//! assert_eq!(root.children.len(), 3);
//!
//! // Read branch lengths from edges
//! for &edge_idx in &root.children {
//!   let edge = &graph.edges[edge_idx];
//!   let child = &graph.nodes[edge.child];
//!   println!("{}: {:?}", child.name.as_deref().unwrap_or(""), edge.data.branch_length);
//! }
//! ```
//!
//! # Write a tree
//!
//! ```
//! use util_newick::{newick_from_string, newick_to_string, NewickWriteOptions, NwkStyle};
//!
//! let graph = newick_from_string("(A[&prob=0.95]:0.1,B:0.2);").unwrap();
//!
//! // Plain: strips all annotations
//! let plain = newick_to_string(&graph, &NewickWriteOptions { style: NwkStyle::Plain, ..Default::default() }).unwrap();
//! assert_eq!(plain, "(A:0.1,B:0.2);");
//!
//! // BEAST: preserves annotations in [&k=v] format
//! let beast = newick_to_string(&graph, &NewickWriteOptions::default()).unwrap();
//! assert!(beast.contains("[&prob="));
//! ```
//!
//! # Read annotations
//!
//! The parser auto-detects annotation dialect per comment:
//!
//! ```
//! use util_newick::{newick_from_string, NewickValue};
//!
//! // BEAST-style: [&key=value] -- typed values
//! let g = newick_from_string("(A[&prob=0.95,fixed=TRUE,hpd={1.0,2.0}],B);").unwrap();
//! let a = g.nodes.iter().find(|n| n.name.as_deref() == Some("A")).unwrap();
//! assert!(matches!(a.node_attrs["prob"], NewickValue::Number(_)));
//! assert!(matches!(a.node_attrs["fixed"], NewickValue::Boolean(true)));
//! assert!(matches!(a.node_attrs["hpd"], NewickValue::Array(_)));
//!
//! // NHX-style: [&&NHX:key=value] -- all values are strings
//! let g = newick_from_string("(A[&&NHX:S=human:B=90],B);").unwrap();
//! let a = g.nodes.iter().find(|n| n.name.as_deref() == Some("A")).unwrap();
//! assert!(matches!(&a.node_attrs["S"], NewickValue::String(s) if s == "human"));
//!
//! // Plain comments [text] without & prefix go to raw_comments
//! let g = newick_from_string("(A[a note],B);").unwrap();
//! let a = g.nodes.iter().find(|n| n.name.as_deref() == Some("A")).unwrap();
//! assert_eq!(a.raw_comments, vec!["[a note]"]);
//! ```
//!
//! # Build a tree programmatically
//!
//! ```
//! use util_newick::*;
//! use std::collections::BTreeMap;
//!
//! let mut graph = NewickGraph::new();
//!
//! let root = graph.add_node(NewickNodeData::new().with_name("root"));
//! graph.root = root;
//!
//! let leaf_a = graph.add_node(NewickNodeData::new().with_name("A"));
//! let leaf_b = graph.add_node(NewickNodeData::new().with_name("B"));
//!
//! graph.add_edge(root, leaf_a, NewickEdgeData::new().with_length(0.1));
//! graph.add_edge(root, leaf_b, NewickEdgeData::new().with_length(0.2));
//!
//! let nwk = newick_to_string(&graph, &NewickWriteOptions::default()).unwrap();
//! assert_eq!(nwk, "(A:0.1,B:0.2)root;");
//! ```
//!
//! # Nexus files
//!
//! ```
//! use util_newick::{nexus_from_string, nexus_to_string, NewickWriteOptions};
//!
//! let nexus = "#NEXUS\n\
//!   Begin Trees;\n\
//!     Translate\n\
//!       1 Homo_sapiens,\n\
//!       2 Pan_troglodytes\n\
//!     ;\n\
//!     Tree primates = (1:0.1,2:0.2);\n\
//!   End;\n";
//!
//! let trees = nexus_from_string(nexus).unwrap();
//! assert_eq!(trees.len(), 1);
//! assert_eq!(trees[0].name, "primates");
//!
//! // TRANSLATE tokens are resolved -- node names are full taxon names
//! let names: Vec<_> = trees[0].graph.nodes.iter()
//!   .filter_map(|n| n.name.as_deref())
//!   .collect();
//! assert!(names.contains(&"Homo_sapiens"));
//!
//! // Round-trip back to Nexus
//! let output = nexus_to_string(&trees, &NewickWriteOptions::default()).unwrap();
//! assert!(output.starts_with("#NEXUS"));
//! ```
//!
//! # eNewick networks
//!
//! Hybrid nodes (`#H1`, `#LGT2`, `##R3` acceptor) are detected in labels and
//! merged into a single graph node with multiple parent edges:
//!
//! ```
//! use util_newick::newick_from_string;
//!
//! let g = newick_from_string("((A,(B)x#H1)c,(x#H1,C)d);").unwrap();
//! let hybrid = g.nodes.iter().find(|n| n.hybrid.is_some()).unwrap();
//! assert_eq!(hybrid.name.as_deref(), Some("x"));
//! let h = hybrid.hybrid.as_ref().unwrap();
//! assert_eq!(h.kind.as_deref(), Some("H"));
//! assert_eq!(h.index, 1);
//!
//! // The hybrid node has two parent edges (from c and d)
//! let parent_edges: Vec<_> = g.edges.iter().filter(|e| e.child == g.nodes.iter().position(|n| n.hybrid.is_some()).unwrap()).collect();
//! assert_eq!(parent_edges.len(), 2);
//! ```
//!
//! # Equality
//!
//! `PartialEq` compares trees order-insensitively -- `(A,B)` equals `(B,A)`.
//! Use `eq_ordered` when child order matters:
//!
//! ```
//! use util_newick::newick_from_string;
//!
//! let g1 = newick_from_string("(A,B,C);").unwrap();
//! let g2 = newick_from_string("(C,A,B);").unwrap();
//! assert_eq!(g1, g2);           // order-insensitive
//! assert!(!g1.eq_ordered(&g2)); // order-sensitive
//! ```
//!
//! # Write styles
//!
//! | Style | Output | Compatible with |
//! |-------|--------|-----------------|
//! | `Plain` | name:length only | all tools |
//! | `Beast` (default) | `[&k=v]` annotations | BEAST, FigTree, MrBayes, IQ-TREE, DendroPy |
//! | `Nhx` | `[&&NHX:k=v]` annotations | Forester, ETE, DendroPy |
//!
//! # I/O with files
//!
//! ```no_run
//! use util_newick::{newick_from_reader, newick_to_writer, NewickWriteOptions};
//! use std::fs::File;
//! use std::io::BufReader;
//!
//! // Read
//! let file = File::open("tree.nwk").unwrap();
//! let graph = newick_from_reader(BufReader::new(file)).unwrap();
//!
//! // Write
//! let mut out = File::create("output.nwk").unwrap();
//! newick_to_writer(&mut out, &graph, &NewickWriteOptions::default()).unwrap();
//! ```

pub mod nexus;
pub mod parse;
pub mod types;
pub mod write;

pub use crate::nexus::{nexus_from_reader, nexus_from_string, nexus_to_string, nexus_to_writer};
pub use crate::parse::{newick_from_reader, newick_from_string};
pub use crate::types::{
  NewickEdgeData, NewickEdgeEntry, NewickGraph, NewickHybrid, NewickNodeData, NewickValue, NewickWriteOptions,
  NexusTree, NwkStyle,
};
pub use crate::write::{
  needs_quoting, newick_to_string, newick_to_writer, write_beast_attrs, write_label, write_nhx_attrs,
};

#[cfg(test)]
mod __tests__;

#[cfg(test)]
mod tests {
  use ctor::ctor;

  #[ctor]
  fn init() {
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
