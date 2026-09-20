use crate::timetree::timetree_state::TimetreeState;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;

#[derive(Clone, Default)]
struct RelaxedClockCoeffs {
  k1: f64,
  k2: f64,
}

pub fn apply_relaxed_clock(
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  params: &[f64],
  one_mutation: f64,
  clock_rate: f64,
  state: &mut TimetreeState,
) -> Result<(), Report> {
  let slack = params.first().copied().unwrap_or(1.0);
  let coupling = params.get(1).copied().unwrap_or(1.0);

  let c = 1.0 / one_mutation;

  let mut coeffs: BTreeMap<GraphNodeKey, RelaxedClockCoeffs> = BTreeMap::new();

  graph.iter_depth_first_postorder_forward(|node| {
    let mut node_coeffs = RelaxedClockCoeffs::default();

    let (opt_len, act_len) = if node.is_root {
      (one_mutation, one_mutation)
    } else if let Some(&edge_key) = node.parent_edge_keys.first() {
      let opt_len = branch_lengths[&edge_key].unwrap_or(0.0);
      let act_len = state
        .edge(edge_key)
        .time_length
        .map_or(opt_len, |time_length| time_length * clock_rate);
      (opt_len, act_len)
    } else {
      (one_mutation, one_mutation)
    };

    let denom = opt_len + one_mutation;
    node_coeffs.k2 = slack + c * act_len * act_len / denom;
    node_coeffs.k1 = -2.0 * (c * act_len * opt_len / denom + slack);

    for (child_key, _edge_key) in &node.child_keys {
      if let Some(child_coeffs) = coeffs.get(child_key) {
        let denom = coupling + child_coeffs.k2;
        if denom.abs() > 1e-10 {
          let ratio = coupling / denom;
          node_coeffs.k2 += coupling * (1.0 - ratio).powi(2) + child_coeffs.k2 * ratio.powi(2);
          node_coeffs.k1 += coupling * (1.0 - ratio) * child_coeffs.k1 / denom
            - coupling * child_coeffs.k1 * child_coeffs.k2 / denom.powi(2)
            + coupling * child_coeffs.k1 / denom;
        }
      }
    }

    coeffs.insert(node.key, node_coeffs);
    Ok(())
  })?;

  let mut gammas: BTreeMap<GraphNodeKey, f64> = BTreeMap::new();

  graph.iter_depth_first_preorder_forward(|node| {
    let node_coeffs = coeffs.get(&node.key).cloned().unwrap_or_default();

    let gamma = if node.is_root {
      if node_coeffs.k2.abs() > 1e-10 {
        (-0.5 * node_coeffs.k1 / node_coeffs.k2).max(0.1)
      } else {
        1.0
      }
    } else {
      let parent_gamma = node
        .parent_keys
        .first()
        .and_then(|(parent_key, _)| gammas.get(parent_key))
        .copied()
        .unwrap_or(1.0);

      let denom = coupling + node_coeffs.k2;
      if denom.abs() > 1e-10 {
        ((coupling * parent_gamma - 0.5 * node_coeffs.k1) / denom).max(0.1)
      } else {
        1.0
      }
    };

    gammas.insert(node.key, gamma);
    Ok(())
  })?;

  for (node_key, gamma) in &gammas {
    if let Some(node) = graph.get_node(*node_key) {
      for (_, edge) in graph.parents_of(node) {
        let edge_key = edge.key();
        state.edge_mut(edge_key).gamma = *gamma;
      }
    }
  }

  Ok(())
}
