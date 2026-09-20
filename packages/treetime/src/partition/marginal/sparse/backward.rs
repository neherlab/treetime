use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::shared::update::MarginalBackward;
use crate::partition::marginal::sparse::message::{combine_messages, propagate_raw, propagate_raw_per_site};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::{
  SparseEdgeBackward, SparseEdgeObs, SparseNodeState, SparseSeqDistribution, VarPos,
};
use eyre::Report;
use maplit::btreemap;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPass, GraphPassBackwardContext, GraphPassNodeOutput};
use treetime_primitives::LogLh;
use treetime_utils::interval::range::range_contains;

pub fn process_backward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_states: &BTreeMap<GraphNodeKey, SparseNodeState>,
) -> Result<MarginalBackward<SparseNodeState, SparseEdgeBackward>, Report> {
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_backward(
    node_states,
    &partition.obs_edges,
    |key| treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass"),
    |context| process_node_backward_indexed(partition, gtr, branch_lengths, &context),
  )?;
  Ok(MarginalBackward {
    node_states: outputs.nodes,
    backward: outputs.edges,
  })
}

#[allow(
  clippy::expect_used,
  clippy::unwrap_used,
  reason = "expect on a value an upstream invariant guarantees is present; unwrap on a value an upstream invariant guarantees is present"
)]
fn process_node_backward_indexed(
  partition: &PartitionMarginalSparse,
  gtr: &GTR,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: &GraphPassBackwardContext<'_, SparseNodeState, SparseEdgeObs, SparseNodeState, SparseEdgeBackward>,
) -> Result<GraphPassNodeOutput<SparseNodeState, SparseEdgeBackward>, Report> {
  let alphabet = &partition.alphabet;
  let length = partition.length;
  let obs = &partition.obs_nodes[&context.key];
  let node = context.input.clone();

  let msg_to_parent = if context.is_leaf {
    let fixed = alphabet
      .determined()
      .map(|state| Ok((state, alphabet.get_profile(state)?.clone())))
      .collect::<Result<_, Report>>()?;
    let variable = obs
      .fitch
      .variable
      .iter()
      .map(|(pos, profile)| {
        let dis = alphabet.construct_profile(profile.chars()).unwrap();
        let state = obs
          .fitch
          .chosen_state
          .get(pos)
          .copied()
          .filter(|state| alphabet.is_canonical(*state))
          .unwrap_or_else(|| profile.get_one());
        (*pos, VarPos { dis, state })
      })
      .collect();
    SparseSeqDistribution {
      fixed_counts: obs.composition.clone(),
      variable,
      variable_indel: BTreeSet::new(),
      fixed,
      log_lh: LogLh::ZERO,
    }
  } else {
    let mut variable_pos = btreemap! {};

    let n_children = context.children.len();
    let mut child_states = vec![btreemap! {}; n_children];
    let mut child_messages = Vec::with_capacity(n_children);

    for (ci, child) in context.children.iter().enumerate() {
      let child_edge_obs = &partition.obs_edges[&child.edge_key];
      let child_backward = child
        .edge
        .expect("Backward child edge message must be published before its parent");
      for mutation in child_edge_obs.fitch_subs() {
        variable_pos.insert(mutation.pos(), mutation.reff());
        child_states[ci].insert(mutation.pos(), mutation.qry());
      }
      for (pos, profile) in &child_backward.msg_from_child.variable {
        variable_pos.entry(*pos).or_insert(profile.state);
      }
      child_messages.push(child_backward.msg_from_child.clone());
    }

    for (ci, child) in context.children.iter().enumerate() {
      let child_obs = &partition.obs_nodes[&child.node_key];
      for (pos, parent_state) in &variable_pos {
        if child_states[ci].contains_key(pos) {
          continue;
        }
        let state = if range_contains(&child_obs.non_char, *pos) {
          if range_contains(&child_obs.gaps, *pos) {
            alphabet.gap()
          } else {
            alphabet.unknown()
          }
        } else {
          *parent_state
        };
        child_states[ci].insert(*pos, state);
      }
    }

    combine_messages(
      &obs.composition,
      &child_messages,
      &variable_pos,
      &child_states,
      alphabet,
      context.is_root.then_some(&gtr.pi),
    )?
  };

  let mut node = node;
  let parent_message = if context.is_root {
    node.profile = msg_to_parent;
    None
  } else {
    let (edge_key, edge_obs) = context.parent_edge.expect("Non-root node must own its parent edge");
    let branch_length = fix_branch_length(length, branch_lengths[&edge_key]);
    let msg_from_child = if gtr.has_site_rates() {
      propagate_raw_per_site(
        gtr,
        branch_length,
        true,
        &msg_to_parent,
        edge_obs.transmission.as_deref(),
      )
    } else {
      propagate_raw(
        &gtr.expQt(branch_length).t().to_owned(),
        &msg_to_parent,
        edge_obs.transmission.as_deref(),
      )
    };
    Some(SparseEdgeBackward {
      msg_to_parent,
      msg_from_child,
    })
  };

  Ok(GraphPassNodeOutput { node, parent_message })
}
