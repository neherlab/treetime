use crate::edge::GraphEdge;
use crate::graph::Graph;
use crate::node::{GraphNode, GraphNodeKey};
use eyre::{Report, eyre};
use itertools::Itertools;

/// Most recent common ancestor of a set of nodes.
///
/// Walks the root-to-node paths of every input node in parallel and returns the
/// last node shared by all of them. With a single input the node itself is its
/// own ancestor. Errors on an empty input set.
pub fn common_ancestor<N, E, D>(graph: &Graph<N, E, D>, node_keys: &[GraphNodeKey]) -> Result<GraphNodeKey, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: Send + Sync,
{
  let paths = node_keys
    .iter()
    .map(|key| {
      graph
        .path_from_root_to_node(*key)
        .map(|path| path.into_iter().map(|node| node.read_arc().key()).collect_vec())
    })
    .try_collect::<_, Vec<_>, _>()?;

  let first_path = paths
    .first()
    .ok_or_else(|| eyre!("Cannot find MRCA of an empty node set"))?;

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
