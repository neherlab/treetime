use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::indexed_pass::{IndexedPassDependencies, IndexedPassSlot};
use crate::partition::marginal::sparse::message::{combine_messages, propagate_raw, propagate_raw_per_site};
use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition, SparseSeqDistribution, VarPos};
use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeSet;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_primitives::LogLh;
use treetime_utils::interval::range::range_contains;

#[allow(clippy::too_many_arguments)]
pub(crate) fn process_node_backward_indexed<N, E>(
  graph: &Graph<N, E, ()>,
  alphabet: &crate::alphabet::alphabet::Alphabet,
  gtr: &crate::gtr::gtr::GTR,
  length: usize,
  dependencies: &IndexedPassDependencies<SparseNodePartition, SparseEdgePartition>,
  slot: &mut IndexedPassSlot<SparseNodePartition, SparseEdgePartition>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let graph_node = graph.get_node(slot.key).expect("Indexed node must exist in graph");
  let graph_node = graph_node.read_arc();
  let msg_to_parent = if graph_node.is_leaf() {
    let fixed = alphabet
      .determined()
      .map(|state| Ok((state, alphabet.get_profile(state)?.clone())))
      .collect::<Result<_, Report>>()?;
    let variable = slot
      .node
      .seq
      .fitch
      .variable
      .iter()
      .map(|(pos, profile)| {
        let dis = alphabet.construct_profile(profile.chars()).unwrap();
        let state = slot
          .node
          .seq
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
      fixed_counts: slot.node.seq.composition.clone(),
      variable,
      variable_indel: BTreeSet::new(),
      fixed,
      log_lh: LogLh::ZERO,
    }
  } else {
    let mut variable_pos = btreemap! {};
    let child_pairs = graph.children_of(&graph_node);
    let mut child_states = vec![btreemap! {}; child_pairs.len()];
    let mut child_messages = Vec::with_capacity(child_pairs.len());
    let mut child_keys = Vec::with_capacity(child_pairs.len());

    for (ci, (child, edge)) in child_pairs.iter().enumerate() {
      let child_key = child.read_arc().key();
      let edge_key = edge.read_arc().key();
      let edge_data = dependencies.edge(edge_key);
      for mutation in edge_data.fitch_subs() {
        variable_pos.insert(mutation.pos(), mutation.reff());
        child_states[ci].insert(mutation.pos(), mutation.qry());
      }
      for (pos, profile) in &edge_data.msg_from_child.variable {
        variable_pos.entry(*pos).or_insert(profile.state);
      }
      child_messages.push(edge_data.msg_from_child.clone());
      child_keys.push(child_key);
    }

    for (ci, child_key) in child_keys.iter().enumerate() {
      let child_data = dependencies.node(*child_key);
      for (pos, parent_state) in &variable_pos {
        if child_states[ci].contains_key(pos) {
          continue;
        }
        let state = if range_contains(&child_data.seq.non_char, *pos) {
          if range_contains(&child_data.seq.gaps, *pos) {
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
      &slot.node.seq.composition,
      &child_messages,
      &variable_pos,
      &child_states,
      alphabet,
      graph_node.is_root().then_some(&gtr.pi),
    )?
  };

  if graph_node.is_root() {
    slot.node.profile = msg_to_parent;
  } else {
    let (edge_key, edge_data) = slot
      .parent_edge
      .as_mut()
      .expect("Non-root node must own its parent edge");
    let branch_length = graph
      .get_edge(*edge_key)
      .expect("Indexed edge must exist in graph")
      .read_arc()
      .payload()
      .read_arc()
      .profile_branch_length()
      .unwrap_or(0.0);
    let branch_length = fix_branch_length(length, branch_length);
    edge_data.msg_from_child = if gtr.has_site_rates() {
      propagate_raw_per_site(
        gtr,
        branch_length,
        true,
        &msg_to_parent,
        edge_data.transmission.as_deref(),
      )
    } else {
      propagate_raw(
        &gtr.expQt(branch_length).t().to_owned(),
        &msg_to_parent,
        edge_data.transmission.as_deref(),
      )
    };
    edge_data.msg_to_parent = msg_to_parent;
  }
  Ok(())
}
