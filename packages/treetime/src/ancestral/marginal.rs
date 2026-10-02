use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;

pub fn branch_lengths_or_zero(branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> BTreeMap<GraphEdgeKey, f64> {
  branch_lengths
    .iter()
    .map(|(key, raw)| (*key, raw.unwrap_or(0.0)))
    .collect()
}
