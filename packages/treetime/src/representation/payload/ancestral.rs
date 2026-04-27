use crate::o;
use crate::representation::partition::traits::PartitionBranchOps;
use eyre::Report;
use itertools::Itertools;
use parking_lot::RwLock;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{Described, GraphNode, Named};
use treetime_io::graphviz::{EdgeToGraphviz, NodeToGraphviz};
use treetime_io::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions, format_weight};

pub type GraphAncestral = Graph<NodeAncestral, EdgeAncestral, ()>;

/// Populate the Newick/Nexus branch-mutations comment on a node payload.
///
/// Implemented by any tree-node payload that carries branch mutation
/// annotations (currently `NodeAncestral` directly, and `NodeTimetree` via
/// its inner `NodeAncestral`). Lets `annotate_branch_mutations()` stay
/// generic over the concrete graph payload used by a given command.
pub trait HasBranchMutations {
  fn set_branch_mutations(&mut self, mutations: Option<String>);
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct NodeAncestral {
  pub name: Option<String>,
  pub desc: Option<String>,
  /// Branch mutations annotation for Newick/Nexus output.
  /// Populated by `annotate_branch_mutations()` from partition data
  /// before tree output. Format: "A55G,T93C" (1-based positions).
  #[serde(skip_serializing_if = "Option::is_none")]
  pub mutations: Option<String>,
}

impl HasBranchMutations for NodeAncestral {
  fn set_branch_mutations(&mut self, mutations: Option<String>) {
    self.mutations = mutations;
  }
}

impl NodeFromNwk for NodeAncestral {
  fn from_nwk(name: Option<impl AsRef<str>>, _: &BTreeMap<String, String>) -> Result<Self, Report> {
    Ok(Self {
      name: name.map(|s| s.as_ref().to_owned()),
      ..NodeAncestral::default()
    })
  }
}

impl NodeToNwk for NodeAncestral {
  fn nwk_name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn nwk_comments(&self) -> BTreeMap<String, String> {
    let mut comments = BTreeMap::new();
    if let Some(mutations) = &self.mutations {
      comments.insert(o!("mutations"), mutations.clone());
    }
    comments
  }
}

impl GraphNode for NodeAncestral {}

impl Named for NodeAncestral {
  fn name(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }

  fn set_name(&mut self, name: Option<impl AsRef<str>>) {
    self.name = name.map(|n| n.as_ref().to_owned());
  }
}

impl Described for NodeAncestral {
  fn desc(&self) -> &Option<String> {
    &self.desc
  }

  fn set_desc(&mut self, desc: Option<String>) {
    self.desc = desc;
  }
}

impl NodeToGraphviz for NodeAncestral {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self.name.as_deref()
  }
}

#[derive(Clone, Default, Debug, Serialize, Deserialize)]
pub struct EdgeAncestral {
  pub branch_length: Option<f64>,
}

impl GraphEdge for EdgeAncestral {}

impl HasBranchLength for EdgeAncestral {
  fn branch_length(&self) -> Option<f64> {
    self.branch_length
  }

  fn set_branch_length(&mut self, weight: Option<f64>) {
    self.branch_length = weight;
  }
}

impl EdgeFromNwk for EdgeAncestral {
  fn from_nwk(branch_length: Option<f64>) -> Result<Self, Report> {
    Ok(Self { branch_length })
  }
}

impl EdgeToNwk for EdgeAncestral {
  fn nwk_weight(&self) -> Option<f64> {
    self.branch_length()
  }
}

impl EdgeToGraphviz for EdgeAncestral {
  fn to_graphviz_label(&self) -> Option<impl AsRef<str>> {
    self
      .branch_length()
      .map(|weight| format_weight(weight, &NwkWriteOptions::default()))
  }

  fn to_graphviz_weight(&self) -> Option<f64> {
    self.branch_length()
  }
}

/// Populate mutation annotations on graph nodes from partition data.
///
/// For each edge, derives current mutations via `PartitionBranchOps::edge_subs()`
/// across all partitions and stores the formatted string on the child node.
/// Mutations are sorted by position and formatted as "A55G,T93C" (1-based).
///
/// Generic over the graph payload: any graph whose node type implements
/// `HasBranchMutations` (currently `NodeAncestral`, and `NodeTimetree` via
/// its inner `NodeAncestral`) can be annotated. The partition slice is
/// `&[Arc<RwLock<P>>]` with `P: PartitionBranchOps + ?Sized`, which accepts
/// both plain `dyn PartitionBranchOps` trait objects and richer trait
/// objects that include `PartitionBranchOps` in their super-trait set
/// (e.g. `dyn PartitionTimetreeAll<N, E>`).
///
/// Call this after reconstruction and before Newick/Nexus output to fill
/// the `mutations` field in `NodeAncestral::nwk_comments()`.
pub fn annotate_branch_mutations<N, E, D, P>(
  graph: &Graph<N, E, D>,
  partitions: &[Arc<RwLock<P>>],
) -> Result<(), Report>
where
  N: GraphNode + HasBranchMutations,
  E: GraphEdge,
  D: Send + Sync,
  P: PartitionBranchOps + ?Sized,
{
  for edge in graph.get_edges() {
    let edge = edge.read_arc();
    let edge_key = edge.key();
    let child_key = edge.target();

    let mut all_subs = Vec::new();
    for partition in partitions {
      let partition = partition.read_arc();
      all_subs.extend(partition.edge_subs(graph, edge_key)?);
    }

    all_subs.sort_by_key(|s| s.pos());

    let mutations = if all_subs.is_empty() {
      None
    } else {
      Some(all_subs.iter().join(","))
    };

    let child = graph.get_node(child_key).expect("Child node must exist");
    child.read_arc().payload().write_arc().set_branch_mutations(mutations);
  }

  Ok(())
}

#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
  use crate::representation::partition::traits::PartitionBranchOps;
  use crate::representation::payload::ancestral::{GraphAncestral, annotate_branch_mutations};
  use crate::representation::payload::sparse::{SparseEdgePartition, SparseNodePartition};
  use crate::seq::mutation::Sub;
  use eyre::Report;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::AsciiChar;

  fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  fn make_test_partition(
    graph: &GraphAncestral,
    length: usize,
    edge_subs: &[(usize, Vec<Sub>)],
  ) -> Result<Arc<RwLock<PartitionMarginalSparse>>, Report> {
    let alphabet = Alphabet::default();
    // Build reference sequence consistent with sub ref chars
    let mut ref_seq: treetime_primitives::Seq = std::iter::repeat_with(|| c(b'A')).take(length).collect();
    for (_, subs) in edge_subs {
      for s in subs {
        if s.pos() < length {
          ref_seq[s.pos()] = s.reff();
        }
      }
    }

    let mut partition = PartitionMarginalSparse {
      index: 0,
      gtr: jc69(JC69Params::default())?,
      alphabet: alphabet.clone(),
      length,
      root_sequence: ref_seq.clone(),
      nodes: btreemap! {},
      edges: btreemap! {},
    };

    for node in graph.get_nodes() {
      let key = node.read_arc().key();
      let mut node_part = SparseNodePartition::empty(&alphabet);
      node_part.seq.sequence = ref_seq.clone();
      partition.nodes.insert(key, node_part);
    }

    let edges = graph.get_edges();
    for (idx, subs) in edge_subs {
      if let Some(edge) = edges.get(*idx) {
        let edge_key = edge.read_arc().key();
        let mut edge_part = SparseEdgePartition::with_fitch_subs(subs.clone());
        edge_part.set_ml_subs(subs.clone());
        partition.edges.insert(edge_key, edge_part);
      }
    }

    Ok(Arc::new(RwLock::new(partition)))
  }

  #[test]
  fn test_annotate_branch_mutations_formats_1_based_positions() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partition = make_test_partition(
      &graph,
      100,
      &[(
        0,
        vec![
          Sub::new(c(b'A'), 0_usize, c(b'T'))?,
          Sub::new(c(b'G'), 5_usize, c(b'C'))?,
        ],
      )],
    )?;
    let branch_ops: Vec<Arc<RwLock<dyn PartitionBranchOps>>> = vec![partition];
    annotate_branch_mutations(&graph, &branch_ops)?;

    let child = graph.get_leaves()[0].read_arc();
    let mutations = child.payload().read_arc().mutations.clone();
    // Sub::Display uses 1-based positions: A1T, G6C
    assert_eq!(mutations, Some("A1T,G6C".to_owned()));
    Ok(())
  }

  #[test]
  fn test_annotate_branch_mutations_empty_partitions() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let branch_ops: Vec<Arc<RwLock<dyn PartitionBranchOps>>> = vec![];
    annotate_branch_mutations(&graph, &branch_ops)?;

    let child = graph.get_leaves()[0].read_arc();
    let mutations = child.payload().read_arc().mutations.clone();
    assert_eq!(mutations, None);
    Ok(())
  }

  #[test]
  fn test_annotate_branch_mutations_no_mutations_on_edge() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partition = make_test_partition(&graph, 100, &[(0, vec![])])?;
    let branch_ops: Vec<Arc<RwLock<dyn PartitionBranchOps>>> = vec![partition];
    annotate_branch_mutations(&graph, &branch_ops)?;

    let child = graph.get_leaves()[0].read_arc();
    let mutations = child.payload().read_arc().mutations.clone();
    assert_eq!(mutations, None);
    Ok(())
  }

  #[test]
  fn test_annotate_branch_mutations_sorts_by_position() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;
    let partition = make_test_partition(
      &graph,
      100,
      &[(
        0,
        vec![
          Sub::new(c(b'C'), 50_usize, c(b'G'))?,
          Sub::new(c(b'A'), 10_usize, c(b'T'))?,
          Sub::new(c(b'G'), 30_usize, c(b'C'))?,
        ],
      )],
    )?;
    let branch_ops: Vec<Arc<RwLock<dyn PartitionBranchOps>>> = vec![partition];
    annotate_branch_mutations(&graph, &branch_ops)?;

    let child = graph.get_leaves()[0].read_arc();
    let mutations = child.payload().read_arc().mutations.clone();
    // Sorted by position (1-based): A11T, G31C, C51G
    assert_eq!(mutations, Some("A11T,G31C,C51G".to_owned()));
    Ok(())
  }

  #[test]
  fn test_annotate_branch_mutations_multi_partition_merge() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str("(A:0.1)root;")?;

    let p1 = make_test_partition(&graph, 100, &[(0, vec![Sub::new(c(b'A'), 5_usize, c(b'T'))?])])?;
    let p2 = make_test_partition(&graph, 100, &[(0, vec![Sub::new(c(b'G'), 20_usize, c(b'C'))?])])?;

    let branch_ops: Vec<Arc<RwLock<dyn PartitionBranchOps>>> = vec![p1, p2];
    annotate_branch_mutations(&graph, &branch_ops)?;

    let child = graph.get_leaves()[0].read_arc();
    let mutations = child.payload().read_arc().mutations.clone();
    // Merged from both partitions, sorted by position
    assert_eq!(mutations, Some("A6T,G21C".to_owned()));
    Ok(())
  }
}
