use crate::graph::Graph;
use crate::node::GraphNodeKey;
use eyre::Report;
use itertools::Itertools;
use treetime_utils::make_report;

pub fn common_ancestor(graph: &Graph, node_keys: &[GraphNodeKey]) -> Result<GraphNodeKey, Report> {
  let paths = node_keys
    .iter()
    .map(|key| graph.path_from_root_to_node(*key))
    .try_collect::<_, Vec<_>, _>()?;

  let first_path = paths
    .first()
    .ok_or_else(|| make_report!("Cannot find MRCA of an empty node set"))?;

  let mut ancestor = first_path[0];
  for (index, candidate) in first_path.iter().copied().enumerate() {
    if paths.iter().all(|path| path.get(index).copied() == Some(candidate)) {
      ancestor = candidate;
    } else {
      break;
    }
  }

  Ok(ancestor)
}
