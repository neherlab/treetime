use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::pair_by_name::pair_by_name;

pub fn pair_trait_rows(
  rows: impl IntoIterator<Item = (String, String)>,
  graph: &Graph,
  names: &NodeNames,
) -> TraitPairing {
  let pairing = pair_by_name(graph.get_leaves().map(|leaf| leaf.key()), names, rows);
  let observed_values = pairing
    .by_node
    .values()
    .chain(pairing.unmatched.iter().map(|(_, value)| value))
    .cloned()
    .collect();
  TraitPairing {
    traits: pairing.by_node,
    observed_values,
    duplicate_names: pairing.duplicate_entry_names,
    unmatched: pairing.unmatched,
  }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TraitPairing {
  pub traits: BTreeMap<GraphNodeKey, String>,
  pub observed_values: BTreeSet<String>,
  pub duplicate_names: Vec<String>,
  pub unmatched: Vec<(String, String)>,
}
