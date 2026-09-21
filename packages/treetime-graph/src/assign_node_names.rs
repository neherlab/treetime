use crate::graph::Graph;
use crate::graph_traverse::GraphNodeForward;
use crate::node::GraphNodeKey;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};

pub fn assign_node_names(
  mut names: BTreeMap<GraphNodeKey, Option<String>>,
  graph: &Graph,
) -> Result<BTreeMap<GraphNodeKey, Option<String>>, Report> {
  let mut result: BTreeMap<GraphNodeKey, Option<String>> = BTreeMap::new();
  let mut used: BTreeSet<String> = BTreeSet::new();
  for node in graph.get_nodes() {
    let key = node.key();
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
