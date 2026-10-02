use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;

pub fn branch_lengths_or_zero(branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> BTreeMap<GraphEdgeKey, f64> {
  branch_lengths
    .iter()
    .map(|(key, raw)| (*key, raw.unwrap_or(0.0)))
    .collect()
}

pub(crate) fn branch_length_or_zero(
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  edge_key: GraphEdgeKey,
) -> f64 {
  branch_lengths[&edge_key].unwrap_or(0.0)
}

#[expect(
  clippy::as_conversions,
  reason = "a sequence length is far below 2^53, so the conversion to f64 is exact"
)]
pub(crate) fn one_mutation(sequence_length: usize) -> f64 {
  1.0 / sequence_length as f64
}
