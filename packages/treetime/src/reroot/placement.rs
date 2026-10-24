use crate::error::input_error;
use crate::reroot::find_root::EdgePoint;
use eyre::Report;
use std::collections::BTreeSet;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::split_fraction::SplitFraction;

const EDGE_END_SNAP_FRACTION: f64 = 1e-5;

pub(crate) fn root_moves(graph: &Graph, point: Option<EdgePoint>, split_edge: bool) -> Result<bool, Report> {
  let Some(point) = point else {
    return Ok(false);
  };
  Ok(RootTarget::on_edge(graph, point, split_edge)? != RootTarget::Node(graph.root_key()?))
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum RootTarget {
  Node(GraphNodeKey),
  Split {
    edge_key: GraphEdgeKey,
    split: SplitFraction,
  },
}

impl RootTarget {
  pub(crate) fn on_edge(graph: &Graph, point: EdgePoint, split_edge: bool) -> Result<Self, Report> {
    let EdgePoint { edge: edge_key, split } = point;
    let (source_key, target_key) = graph.edge_endpoints(edge_key)?;
    let x = split.value();
    let at_target = 1.0 - x < EDGE_END_SNAP_FRACTION || (!split_edge && x >= 0.5);
    Ok(if x < EDGE_END_SNAP_FRACTION {
      Self::Node(source_key)
    } else if at_target && graph.is_leaf(target_key) {
      if split_edge {
        Self::Split {
          edge_key,
          split: SplitFraction::CHILD_END,
        }
      } else {
        Self::Node(source_key)
      }
    } else if at_target {
      Self::Node(target_key)
    } else if split_edge {
      Self::Split { edge_key, split }
    } else {
      Self::Node(source_key)
    })
  }
}

pub(crate) fn leaf_keys(graph: &Graph) -> BTreeSet<GraphNodeKey> {
  graph.get_leaves().map(|leaf| leaf.key()).collect()
}

pub(crate) fn require_dated_new_leaves(
  graph: &Graph,
  leaves_before: &BTreeSet<GraphNodeKey>,
  is_dated: impl Fn(GraphNodeKey) -> bool,
  names: &NodeNames,
) -> Result<(), Report> {
  let undated = graph
    .get_leaves()
    .map(|leaf| leaf.key())
    .find(|key| !leaves_before.contains(key) && !is_dated(*key));
  match undated {
    None => Ok(()),
    Some(key) => Err(input_error(format!(
      "Rerooting turned the internal node '{}' into a leaf, but it has no date and no observed data. \
       Remove the node from the input tree or pass --keep-root.",
      &names[key]
    ))),
  }
}
