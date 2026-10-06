use crate::graph::Graph;
use crate::graph_traverse::GraphNodeForward;
use crate::node::GraphNodeKey;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};

pub fn assign_node_names(
  mut names: BTreeMap<GraphNodeKey, Option<String>>,
  graph: &Graph,
) -> Result<AssignedNodeNames, Report> {
  let mut result: BTreeMap<GraphNodeKey, Option<String>> = BTreeMap::new();
  let mut used: BTreeSet<String> = BTreeSet::new();
  let mut duplicate_names: BTreeSet<String> = BTreeSet::new();
  for node in graph.get_nodes() {
    let key = node.key();
    let name = names.remove(&key).flatten().filter(|name| !name.is_empty());
    if let Some(name) = &name
      && !used.insert(name.clone())
    {
      duplicate_names.insert(name.clone());
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

  Ok(AssignedNodeNames {
    names: result,
    duplicate_names: duplicate_names.into_iter().collect(),
  })
}

pub struct AssignedNodeNames {
  pub names: BTreeMap<GraphNodeKey, Option<String>>,
  pub duplicate_names: Vec<String>,
}

pub fn restrict_node_names(
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  graph: &Graph,
) -> BTreeMap<GraphNodeKey, Option<String>> {
  graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      (key, names.get(&key).cloned().flatten())
    })
    .collect()
}

pub fn node_name_or_key(key: GraphNodeKey, name: Option<&str>) -> String {
  name.map_or_else(|| format!("node_{}", key.as_usize()), str::to_owned)
}
