use crate::graph::Graph;
use crate::key_map::NodeMap;
use crate::node::GraphNodeKey;
use crate::node_names::NodeNames;
use crate::schedule::GraphSchedule;
use eyre::Report;
use indexmap::{IndexMap, IndexSet};
use itertools::{Itertools, izip};
use std::collections::BTreeMap;

pub fn clade_keys<'a>(nodes: impl IntoIterator<Item = CladeNode<'a>>, ids: &mut TipIds<'a>) -> Vec<Vec<usize>> {
  let nodes = nodes.into_iter().collect_vec();
  let mut keys = nodes
    .iter()
    .map(|node| node.tip.map_or_else(Vec::new, |name| vec![ids.id(name)]))
    .collect_vec();
  for (index, node) in nodes.iter().enumerate().rev() {
    keys[index].sort_unstable();
    if let Some(parent) = node.parent {
      assert!(
        parent < index,
        "node {index} has parent {parent}, which does not precede it in preorder"
      );
      let (head, tail) = keys.split_at_mut(index);
      head[parent].extend_from_slice(&tail[0]);
    }
  }
  keys
}

pub fn graph_clade_keys<'a>(
  graph: &Graph,
  names: &'a NodeNames,
  ids: &mut TipIds<'a>,
) -> Result<IndexMap<GraphNodeKey, Vec<usize>>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let preorder = schedule.preorder();
  let positions: NodeMap<usize> = preorder.iter().enumerate().map(|(index, &key)| (key, index)).collect();
  let nodes = preorder.iter().map(|&key| CladeNode {
    parent: schedule.parents(key).first().map(|(parent, _)| positions[parent]),
    tip: schedule.children(key).is_empty().then(|| &names[key]),
  });
  let keys = clade_keys(nodes, ids);
  Ok(izip!(preorder.iter().copied(), keys).collect())
}

pub fn first_node_by_key<'k, N>(keys: impl IntoIterator<Item = (N, &'k [usize])>) -> BTreeMap<&'k [usize], N> {
  let mut first = BTreeMap::new();
  for (node, key) in keys {
    first.entry(key).or_insert(node);
  }
  first
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct CladeNode<'a> {
  pub parent: Option<usize>,
  pub tip: Option<&'a str>,
}

#[derive(Clone, Debug, Default)]
pub struct TipIds<'a>(IndexSet<&'a str>);

impl<'a> TipIds<'a> {
  pub fn id(&mut self, name: &'a str) -> usize {
    self.0.insert_full(name).0
  }

  pub fn name(&self, id: usize) -> &'a str {
    self.0[id]
  }
}
