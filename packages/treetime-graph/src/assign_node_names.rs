use crate::edge::GraphEdge;
use crate::graph::Graph;
use crate::graph_traverse::GraphNodeForward;
use crate::node::{GraphNode, GraphNodeKey};
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};

/// Assign synthetic `NODE_{counter:07}` names to unnamed nodes in DFS-preorder and return the
/// completed node-keyed name map.
///
/// The input `names` supplies each node's known name: the parsed name for named nodes and the name a
/// previous call assigned, with `None` (or an empty name) where a node still needs one. The result is
/// keyed by exactly the graph's current nodes: entries for nodes that a topology change removed are
/// dropped, and nodes a topology change introduced (absent from the input map) are named here. Each
/// unnamed node takes the next free `NODE_xxxxx`, skipping any name already held by a current node, so
/// the numbering is deterministic and stable across re-runs after a topology change. Names are a
/// threaded value rather than node payload state.
pub fn assign_node_names<N: GraphNode, E: GraphEdge, D: Sync + Send>(
  mut names: BTreeMap<GraphNodeKey, Option<String>>,
  graph: &Graph<N, E, D>,
) -> Result<BTreeMap<GraphNodeKey, Option<String>>, Report> {
  let mut result: BTreeMap<GraphNodeKey, Option<String>> = BTreeMap::new();
  let mut used: BTreeSet<String> = BTreeSet::new();
  for node in graph.get_nodes() {
    let key = node.read_arc().key();
    let name = names.remove(&key).flatten().filter(|name| !name.is_empty());
    if let Some(name) = &name {
      used.insert(name.clone());
    }
    result.insert(key, name);
  }

  let mut internal_node_counter = 0;

  graph.iter_depth_first_preorder_forward(|GraphNodeForward { key, .. }| {
    if result.get(&key).is_none_or(Option::is_none) {
      let mut name = format!("NODE_{internal_node_counter:07}");
      while used.contains(&name) {
        internal_node_counter += 1;
        name = format!("NODE_{internal_node_counter:07}");
      }
      used.insert(name.clone());
      result.insert(key, Some(name));
      internal_node_counter += 1;
    }

    Ok(())
  })?;

  Ok(result)
}
