use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::indexed_pass::{IndexedPass, IndexedPassSlot};
use crate::partition::marginal_core::{
  forward_log_lh_add_normalization, forward_log_lh_remove_child, normalize_1d_inplace,
};
use crate::partition::marginal_helpers::{combine_messages, propagate_raw, propagate_raw_per_site};
use crate::partition::marginal_sparse::{PartitionMarginalSparse, reconstruct_map_seq};
use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition, SparseSeqDistribution, VarPos};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use parking_lot::RwLock;
use rayon::prelude::*;
use std::collections::BTreeSet;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_primitives::AsciiChar;
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::interval::range::range_contains;

pub fn process_backward_indexed<N, E>(
  partition: &mut PartitionMarginalSparse,
  graph: &Graph<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let alphabet = partition.alphabet.clone();
  let gtr = partition.gtr.clone();
  let length = partition.length;
  let (nodes, edges) = (&mut partition.nodes, &mut partition.edges);
  let mut pass = IndexedPass::new(graph, nodes, edges, |key| {
    treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass")
  })?;
  let result = pass.try_for_each_backward_frontier(|node_indices, edge_indices, edges, _, completed, frontier| {
    frontier.par_iter_mut().try_for_each(|slot| {
      process_node_backward_indexed(
        graph,
        &alphabet,
        &gtr,
        length,
        node_indices,
        edge_indices,
        edges,
        completed,
        slot,
      )
    })
  });
  let (nodes, edges) = pass.into_maps()?;
  partition.nodes = nodes;
  partition.edges = edges;
  result
}

pub fn process_forward_indexed<N, E>(
  partition: &mut PartitionMarginalSparse,
  graph: &Graph<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let alphabet = partition.alphabet.clone();
  let gtr = partition.gtr.clone();
  let length = partition.length;
  let root_sequence = partition.root_sequence.clone();
  let (nodes, edges) = (&mut partition.nodes, &mut partition.edges);
  let mut pass = IndexedPass::new(graph, nodes, edges, |key| {
    treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass")
  })?;
  let result = pass.try_for_each_forward_frontier(
    |node_indices, edge_indices, edges, future, completed_start, frontier, completed| {
      frontier.par_iter_mut().try_for_each(|slot| {
        process_node_forward_indexed(
          graph,
          &alphabet,
          &gtr,
          length,
          &root_sequence,
          node_indices,
          edge_indices,
          edges,
          future,
          completed_start,
          completed,
          slot,
        )
      })
    },
  );
  let (nodes, edges) = pass.into_maps()?;
  partition.nodes = nodes;
  partition.edges = edges;
  result
}

#[allow(clippy::too_many_arguments)]
fn process_node_forward_indexed<N, E>(
  graph: &Graph<N, E, ()>,
  alphabet: &crate::alphabet::alphabet::Alphabet,
  gtr: &crate::gtr::gtr::GTR,
  length: usize,
  root_sequence: &treetime_primitives::Seq,
  node_indices: &[Option<usize>],
  edge_indices: &[Option<usize>],
  edges: &[RwLock<Option<(GraphEdgeKey, SparseEdgePartition)>>],
  future: &[IndexedPassSlot<SparseNodePartition>],
  completed_start: usize,
  completed: &[IndexedPassSlot<SparseNodePartition>],
  slot: &mut IndexedPassSlot<SparseNodePartition>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  if let Some(parent_key) = slot.parent_key {
    let parent_index = node_indices[parent_key.as_usize()].expect("Indexed parent must have a slot");
    let parent = &completed[parent_index - completed_start].node;
    let slot_index = node_indices[slot.key.as_usize()].expect("Indexed node must have a slot");
    let mut edge = edges[slot_index].write();
    let (edge_key, edge_data) = edge.as_mut().expect("Non-root node must own its parent edge");

    let mut variable_pos = btreemap! {};
    let mut parent_state = btreemap! {};
    let mut child_state = btreemap! {};
    for mutation in edge_data.fitch_subs() {
      let current_state = mutation.qry();
      variable_pos.insert(mutation.pos(), current_state);
      parent_state.entry(mutation.pos()).or_insert_with(|| mutation.reff());
      child_state.entry(mutation.pos()).or_insert(current_state);
    }
    for (pos, profile) in &edge_data.msg_to_child.variable {
      if !range_contains(&slot.node.seq.non_char, *pos) {
        variable_pos.entry(*pos).or_insert(profile.state);
        parent_state.entry(*pos).or_insert(profile.state);
      }
    }
    for (pos, profile) in &edge_data.msg_to_parent.variable {
      variable_pos.entry(*pos).or_insert(profile.state);
      child_state.entry(*pos).or_insert(profile.state);
    }

    let branch_length = graph
      .get_edge(*edge_key)
      .expect("Indexed edge must exist in graph")
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap_or(0.0);
    let branch_length = fix_branch_length(length, branch_length);
    let msg_from_parent = if gtr.has_site_rates() {
      propagate_raw_per_site(
        gtr,
        branch_length,
        false,
        &edge_data.msg_to_child,
        edge_data.transmission.as_deref(),
      )
    } else {
      propagate_raw(
        &gtr.expQt(branch_length),
        &edge_data.msg_to_child,
        edge_data.transmission.as_deref(),
      )
    };
    let profile = combine_messages(
      &slot.node.seq.composition,
      &[msg_from_parent, edge_data.msg_to_parent.clone()],
      &variable_pos,
      &[parent_state, child_state],
      alphabet,
      None,
    )?;
    slot.node.profile = profile;

    if !parent.seq.sequence.is_empty() {
      slot.node.seq.sequence = reconstruct_map_seq(&parent.seq.sequence, Some(edge_data), &slot.node, alphabet);
    }
    edge_data.set_ml_subs(compute_ml_subs_for_nodes(alphabet, parent, &slot.node, edge_data)?);
  } else if slot.node.seq.sequence.is_empty() {
    slot.node.seq.sequence = root_sequence.clone();
  }

  let graph_node = graph.get_node(slot.key).expect("Indexed node must exist in graph");
  let graph_node = graph_node.read_arc();
  for (child, child_edge) in graph.children_of(&graph_node) {
    let child_key = child.read_arc().key();
    let child_index = node_indices[child_key.as_usize()].expect("Indexed child must have a slot");
    let child_node = &future[child_index].node;
    let child_edge_key = child_edge.read_arc().key();
    let child_edge_index = edge_indices[child_edge_key.as_usize()].expect("Indexed child edge must have a slot");
    let mut child_edge = edges[child_edge_index].write();
    let (_, child_edge) = child_edge.as_mut().expect("Indexed child edge slot must be populated");

    let mut seq_dis = SparseSeqDistribution {
      variable: btreemap! {},
      variable_indel: BTreeSet::new(),
      fixed: btreemap! {},
      fixed_counts: slot.node.seq.composition.clone(),
      log_lh: forward_log_lh_remove_child(slot.node.profile.log_lh, child_edge.msg_from_child.log_lh),
    };
    let child_dis = child_edge.msg_from_child.clone();
    let mut parent_states = btreemap! {};
    let mut child_states = btreemap! {};
    for mutation in child_edge.fitch_subs() {
      child_states.insert(mutation.pos(), mutation.qry());
      parent_states.insert(mutation.pos(), mutation.reff());
    }
    for (pos, profile) in &slot.node.profile.variable {
      if !range_contains(&child_node.seq.non_char, *pos) {
        child_states.entry(*pos).or_insert(profile.state);
        parent_states.entry(*pos).or_insert(profile.state);
      }
    }
    for (pos, profile) in &child_dis.variable {
      if !range_contains(&child_node.seq.non_char, *pos) {
        child_states.entry(*pos).or_insert(profile.state);
        parent_states.entry(*pos).or_insert(profile.state);
      }
    }

    let mut delta_ll = 0.0;
    for (pos, parent_state) in parent_states {
      let divisor = child_dis
        .variable
        .get(&pos)
        .map_or(&child_dis.fixed[&child_states[&pos]], |distribution| &distribution.dis);
      let numerator = slot
        .node
        .profile
        .variable
        .get(&pos)
        .map_or(&slot.node.profile.fixed[&parent_state], |distribution| {
          &distribution.dis
        });
      let safe_divisor = divisor.mapv(|value| value.max(f64::MIN_POSITIVE));
      let mut dis = numerator / &safe_divisor;
      let normalization = normalize_1d_inplace(&mut dis, 1.0);
      delta_ll = forward_log_lh_add_normalization(delta_ll, normalization);
      seq_dis.variable.insert(
        pos,
        VarPos {
          dis,
          state: parent_state,
        },
      );
      seq_dis.fixed_counts.adjust_count(parent_state, -1);
    }
    for (state, profile) in &slot.node.profile.fixed {
      let safe_fixed = child_dis.fixed[state].mapv(|value| value.max(f64::MIN_POSITIVE));
      let mut dis = profile / &safe_fixed;
      let weight = seq_dis.fixed_counts.get(*state).unwrap() as f64;
      let normalization = normalize_1d_inplace(&mut dis, weight);
      delta_ll = forward_log_lh_add_normalization(delta_ll, normalization);
      seq_dis.fixed.insert(*state, dis);
    }
    seq_dis.log_lh = forward_log_lh_add_normalization(seq_dis.log_lh, delta_ll);
    child_edge.msg_to_child = seq_dis;
  }
  Ok(())
}

fn compute_ml_subs_for_nodes(
  alphabet: &crate::alphabet::alphabet::Alphabet,
  parent: &SparseNodePartition,
  child: &SparseNodePartition,
  edge: &SparseEdgePartition,
) -> Result<Vec<Sub>, Report> {
  let positions = edge
    .fitch_subs()
    .iter()
    .map(Sub::pos)
    .chain(parent.profile.variable.keys().copied())
    .chain(child.profile.variable.keys().copied())
    .sorted()
    .dedup();
  positions
    .filter_map(|pos| {
      let parent_state = resolve_map_state(parent, pos, alphabet);
      let child_state = resolve_map_state(child, pos, alphabet);
      (parent_state != child_state && alphabet.is_canonical(parent_state) && alphabet.is_canonical(child_state))
        .then(|| Sub::new(parent_state, pos, child_state))
    })
    .collect()
}

#[allow(clippy::too_many_arguments)]
fn process_node_backward_indexed<N, E>(
  graph: &Graph<N, E, ()>,
  alphabet: &crate::alphabet::alphabet::Alphabet,
  gtr: &crate::gtr::gtr::GTR,
  length: usize,
  node_indices: &[Option<usize>],
  edge_indices: &[Option<usize>],
  edges: &[RwLock<Option<(GraphEdgeKey, SparseEdgePartition)>>],
  completed: &[IndexedPassSlot<SparseNodePartition>],
  slot: &mut IndexedPassSlot<SparseNodePartition>,
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
      log_lh: 0.0,
    }
  } else {
    let mut variable_pos = btreemap! {};
    let child_pairs = graph.children_of(&graph_node);
    let mut child_states = vec![btreemap! {}; child_pairs.len()];
    let mut child_messages = Vec::with_capacity(child_pairs.len());
    let mut child_edge_guards = Vec::with_capacity(child_pairs.len());

    for (ci, (child, edge)) in child_pairs.iter().enumerate() {
      let child_key = child.read_arc().key();
      let edge_key = edge.read_arc().key();
      let edge_index = edge_indices[edge_key.as_usize()].expect("Indexed child edge must have a slot");
      let edge = edges[edge_index].read();
      let edge_data = &edge.as_ref().expect("Indexed child edge slot must be populated").1;
      for mutation in edge_data.fitch_subs() {
        variable_pos.insert(mutation.pos(), mutation.reff());
        child_states[ci].insert(mutation.pos(), mutation.qry());
      }
      for (pos, profile) in &edge_data.msg_from_child.variable {
        variable_pos.entry(*pos).or_insert(profile.state);
      }
      child_messages.push(edge_data.msg_from_child.clone());
      child_edge_guards.push((child_key, edge));
    }

    for (ci, (child_key, _)) in child_edge_guards.iter().enumerate() {
      let child_index = node_indices[child_key.as_usize()].expect("Indexed child must have a slot");
      let child_data = &completed[child_index].node;
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
    let slot_index = node_indices[slot.key.as_usize()].expect("Indexed node must have a slot");
    let mut edge = edges[slot_index].write();
    let (edge_key, edge_data) = edge.as_mut().expect("Non-root node must own its parent edge");
    let branch_length = graph
      .get_edge(*edge_key)
      .expect("Indexed edge must exist in graph")
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
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

fn resolve_map_state(
  node: &SparseNodePartition,
  pos: usize,
  alphabet: &crate::alphabet::alphabet::Alphabet,
) -> AsciiChar {
  if let Some(var) = node.profile.variable.get(&pos) {
    alphabet.char(argmax_first(&var.dis.view()).unwrap_or(0))
  } else {
    node.seq.sequence.get(pos).copied().unwrap_or_else(|| alphabet.char(0))
  }
}
