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
/// Call this after reconstruction and before Newick/Nexus output to fill
/// the `mutations` field in `NodeAncestral::nwk_comments()`.
pub fn annotate_branch_mutations(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<dyn PartitionBranchOps>>],
) -> Result<(), Report> {
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
    child.read_arc().payload().write_arc().mutations = mutations;
  }

  Ok(())
}
