use crate::partition::timetree::GraphTimetree;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::{HasBranchLength, TimeLength};
use treetime_graph::node::GraphNodeKey;

/// Relaxed clock penalty coefficients computed during postorder pass.
#[derive(Clone, Default)]
struct RelaxedClockCoeffs {
  k1: f64,
  k2: f64,
}

/// Apply relaxed molecular clock allowing branch-specific rate variation.
///
/// Allows mutation rate to vary across the tree (relaxed molecular clock).
/// Penalizes rate changes from branch to branch using autocorrelated model.
///
/// Parameters in `params`:
/// - `params[0]`: slack - penalty for deviation from mean rate (default 1.0)
/// - `params[1]`: coupling - penalty for rate differences between parent and child (default 1.0)
///
/// `clock_rate` converts edge time lengths into substitutions per site so that
/// the actual and sequence-optimal branch lengths use the same units.
///
/// The algorithm:
/// 1. Postorder pass: compute quadratic penalty coefficients (k1, k2) for each node
/// 2. Preorder pass: compute optimal gamma (rate multiplier) for each branch
pub fn apply_relaxed_clock(
  graph: &GraphTimetree,
  params: &[f64],
  one_mutation: f64,
  clock_rate: f64,
) -> Result<(), Report> {
  let slack = params.first().copied().unwrap_or(1.0);
  let coupling = params.get(1).copied().unwrap_or(1.0);

  // Coefficient for converting between branch length and mutation count
  let c = 1.0 / one_mutation;

  // Store coefficients per node during postorder
  let mut coeffs: BTreeMap<GraphNodeKey, RelaxedClockCoeffs> = BTreeMap::new();

  // Postorder pass: compute k1, k2 coefficients from leaves to root
  graph.iter_depth_first_postorder_forward(|node| {
    let mut node_coeffs = RelaxedClockCoeffs::default();

    // Compute branch penalty coefficients
    // For root: use one_mutation as both opt_len and act_len (v0 lines 1096-1105)
    // For non-root: use actual edge lengths
    let (opt_len, act_len) = if node.is_root {
      (one_mutation, one_mutation)
    } else if let Some(parent_edge) = node.parent_edges.first() {
      let opt_len = parent_edge.branch_length().unwrap_or(0.0);
      let act_len = parent_edge
        .time_length()
        .map_or(opt_len, |time_length| time_length * clock_rate);
      (opt_len, act_len)
    } else {
      (one_mutation, one_mutation)
    };

    // Contact term: stiffness * (gamma * act_len - opt_len)^2 + slack * (gamma - 1)^2
    // Expanding: (slack + c * act_len^2 / (opt_len + one_mutation)) * gamma^2
    //          - 2 * (c * act_len * opt_len / (opt_len + one_mutation) + slack) * gamma + C
    // = k2 * gamma^2 + k1 * gamma + C
    let denom = opt_len + one_mutation;
    node_coeffs.k2 = slack + c * act_len * act_len / denom;
    node_coeffs.k1 = -2.0 * (c * act_len * opt_len / denom + slack);

    // Coupling term: sum over children of coupling * (gamma - gamma_child)^2 + Cost_child(gamma_child | gamma)
    // Given gamma, optimal gamma_child = (coupling * gamma - 0.5 * child.k1) / (coupling + child.k2)
    // Substituting back yields additional terms for k1, k2
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
      // Root: gamma = max(0.1, -0.5 * k1 / k2)
      if node_coeffs.k2.abs() > 1e-10 {
        (-0.5 * node_coeffs.k1 / node_coeffs.k2).max(0.1)
      } else {
        1.0
      }
    } else {
      // Non-root: gamma = max(0.1, (coupling * parent_gamma - 0.5 * k1) / (coupling + k2))
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

  // Store gammas in edges: each edge stores the gamma of its child node
  for (node_key, gamma) in &gammas {
    if let Some(node) = graph.get_node(*node_key) {
      let node = node.read_arc();
      for (_, edge) in graph.parents_of(&node) {
        edge.write_arc().payload().write_arc().gamma = *gamma;
      }
    }
  }

  Ok(())
}
