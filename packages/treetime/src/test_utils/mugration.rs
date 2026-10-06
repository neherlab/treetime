use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pair_by_name::pair_by_name;

pub(crate) struct TraitsByNode {
  pub(crate) traits: BTreeMap<GraphNodeKey, String>,
  pub(crate) observed_values: BTreeSet<String>,
}

pub(crate) fn traits_by_node(
  rows: impl IntoIterator<Item = (String, String)>,
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
) -> TraitsByNode {
  let pairing = pair_by_name(graph.get_leaves().map(|leaf| leaf.key()), names, rows);
  let observed_values = pairing
    .by_node
    .values()
    .chain(pairing.unmatched.iter().map(|(_, value)| value))
    .cloned()
    .collect();
  TraitsByNode {
    traits: pairing.by_node,
    observed_values,
  }
}
