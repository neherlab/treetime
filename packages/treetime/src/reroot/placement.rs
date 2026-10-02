use crate::node_label::node_label;
use approx::ulps_eq;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::reroot::StemRemovalInfo;
use treetime_utils::make_error;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum RootTarget {
  Node(GraphNodeKey),
  Split { edge_key: GraphEdgeKey, split: f64 },
}

impl RootTarget {
  pub(crate) fn on_edge(graph: &Graph, edge_key: GraphEdgeKey, split: f64, split_edge: bool) -> Result<Self, Report> {
    let (source_key, target_key) = graph.edge_endpoints(edge_key)?;
    Ok(if ulps_eq!(split, 0.0, max_ulps = 5) {
      Self::Node(source_key)
    } else if ulps_eq!(split, 1.0, max_ulps = 5) {
      Self::Node(target_key)
    } else if split_edge {
      Self::Split { edge_key, split }
    } else if split < 0.5 {
      Self::Node(source_key)
    } else {
      Self::Node(target_key)
    })
  }

  #[must_use]
  pub(crate) fn after_stem_removal(self, stem: &StemRemovalInfo) -> Self {
    match self {
      Self::Split { edge_key, .. } if edge_key == stem.removed_edge_key => Self::Node(stem.new_root_key),
      target => target,
    }
  }
}

pub(crate) fn leaf_keys(graph: &Graph) -> BTreeSet<GraphNodeKey> {
  graph.get_leaves().map(|leaf| leaf.key()).collect()
}

pub(crate) fn require_dated_new_leaves(
  graph: &Graph,
  leaves_before: &BTreeSet<GraphNodeKey>,
  is_dated: impl Fn(GraphNodeKey) -> bool,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> Result<(), Report> {
  let undated = graph
    .get_leaves()
    .map(|leaf| leaf.key())
    .find(|key| !leaves_before.contains(key) && !is_dated(*key));
  match undated {
    None => Ok(()),
    Some(key) => make_error!(
      "Rerooting turned the internal node '{}' into a leaf, but it has no date and no observed data. \
       Remove the node from the input tree or pass --keep-root.",
      node_label(names, key)
    ),
  }
}
