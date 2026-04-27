use crate::hacks::fix_branch_length::fix_branch_length;
use crate::representation::partition::marginal_helpers::{combine_messages, propagate_raw, propagate_raw_per_site};
use crate::representation::partition::marginal_sparse::{PartitionMarginalSparse, reconstruct_map_seq};
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, SparseNodePartition, VarPos};
use crate::seq::mutation::Sub;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use ndarray::Array1;
use std::collections::BTreeMap;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_primitives::AsciiChar;
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range::range_contains;

fn node_reference_state(partition: &PartitionMarginalSparse, node_key: GraphNodeKey, pos: usize) -> Option<AsciiChar> {
  partition
    .nodes
    .get(&node_key)
    .and_then(|node_data| node_data.seq.fitch.chosen_state.get(&pos).copied())
}

fn node_reference_state_or(
  partition: &PartitionMarginalSparse,
  node_key: GraphNodeKey,
  pos: usize,
  fallback: AsciiChar,
) -> AsciiChar {
  node_reference_state(partition, node_key, pos)
    .filter(|state| partition.alphabet.is_canonical(*state))
    .unwrap_or(fallback)
}

fn reconstruct_map_sequence<N, E>(partition: &mut PartitionMarginalSparse, node: &GraphNodeForward<N, E, ()>)
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let (base_seq, edge) = if node.is_root {
    (&partition.root_sequence, None)
  } else {
    let Some((parent_key, edge_key)) = get_exactly_one(&node.parent_keys).ok() else {
      return;
    };
    let Some(parent_data) = partition.nodes.get(parent_key) else {
      return;
    };
    if parent_data.seq.sequence.is_empty() {
      return;
    }
    let Some(edge_data) = partition.edges.get(edge_key) else {
      return;
    };
    (&parent_data.seq.sequence, Some(edge_data))
  };

  let Some(node_data) = partition.nodes.get(&node.key) else {
    return;
  };

  let seq = reconstruct_map_seq(base_seq, edge, node_data, &partition.alphabet);
  partition.nodes.get_mut(&node.key).unwrap().seq.sequence = seq;
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

fn compute_ml_subs_for_edge(
  partition: &PartitionMarginalSparse,
  parent_key: GraphNodeKey,
  child_key: GraphNodeKey,
  edge_key: GraphEdgeKey,
) -> Result<Vec<Sub>, Report> {
  let alphabet = &partition.alphabet;
  let Some(edge) = partition.edges.get(&edge_key) else {
    return Ok(vec![]);
  };
  let Some(parent_node) = partition.nodes.get(&parent_key) else {
    return Ok(vec![]);
  };
  let Some(child_node) = partition.nodes.get(&child_key) else {
    return Ok(vec![]);
  };

  // Candidate positions: union of Fitch subs, parent variable sites, child variable sites
  let positions: Vec<usize> = edge
    .fitch_subs()
    .iter()
    .map(Sub::pos)
    .chain(parent_node.profile.variable.keys().copied())
    .chain(child_node.profile.variable.keys().copied())
    .sorted()
    .dedup()
    .collect();

  let mut subs = Vec::new();

  for pos in positions {
    let parent_state = resolve_map_state(parent_node, pos, alphabet);
    let child_state = resolve_map_state(child_node, pos, alphabet);

    if parent_state == child_state {
      continue;
    }
    if !alphabet.is_canonical(parent_state) || !alphabet.is_canonical(child_state) {
      continue;
    }

    subs.push(Sub::new(parent_state, pos, child_state)?);
  }

  Ok(subs)
}

pub fn process_node_backward<N, E>(
  partition: &mut PartitionMarginalSparse,
  node: &GraphNodeBackward<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let length = partition.length;

  let msg_to_parent = if node.is_leaf {
    let alphabet = &partition.alphabet;
    let fixed = alphabet
      .determined()
      .map(|state| Ok((state, alphabet.get_profile(state)?.clone())))
      .collect::<Result<_, Report>>()?;

    let node_data = &partition.nodes[&node.key];

    // convert the parsimony variable states to sparseSeqDis format
    let variable = node_data
      .seq
      .fitch
      .variable
      .iter()
      .map(|(pos, p)| {
        let dis = alphabet.construct_profile(p.chars()).unwrap();
        let state = node_reference_state_or(partition, node.key, *pos, p.get_one());
        (*pos, VarPos { dis, state })
      })
      .collect();

    MarginalSparseSeqDistribution {
      fixed_counts: node_data.seq.composition.clone(),
      variable,
      variable_indel: btreemap! {},
      fixed,
      log_lh: 0.0,
    }
  } else {
    let mut variable_pos = btreemap! {};
    let mut child_states = vec![];
    let mut child_messages: Vec<MarginalSparseSeqDistribution> = vec![];
    for (ci, (_child_key, edge_key)) in node.child_keys.iter().enumerate() {
      child_states.push(btreemap! {});
      let edge_data = &partition.edges[edge_key];
      for m in edge_data.fitch_subs() {
        variable_pos.insert(m.pos(), m.reff());
        child_states[ci].insert(m.pos(), m.qry());
      }
      for (pos, p) in &edge_data.msg_from_child.variable {
        // if position in variable in child and we have not yet encountered it (no sub at this position), then state is the same as in the child.
        variable_pos.entry(*pos).or_insert(p.state);
      }
      child_messages.push(edge_data.msg_from_child.clone());
    }

    // Fill in child states for variable positions without edge substitutions.
    let alphabet = &partition.alphabet;
    for (ci, (child_key, _)) in node.child_keys.iter().enumerate() {
      let states = &mut child_states[ci];
      let child_data = &partition.nodes[child_key];
      for (pos, parent_state) in &variable_pos {
        // we already know the state
        if states.contains_key(pos) {
          continue;
        }

        // not knowing the state implies no substitution at this position. can use the child state as the parent state.
        if range_contains(&child_data.seq.non_char, *pos) {
          let state = if range_contains(&child_data.seq.gaps, *pos) {
            alphabet.gap()
          } else {
            alphabet.unknown()
          };
          states.insert(*pos, state);
        } else {
          states.insert(*pos, *parent_state);
        }
      }
    }

    let composition = partition.nodes[&node.key].seq.composition.clone();
    combine_messages(
      &composition,
      &child_messages,
      &variable_pos,
      &child_states,
      &partition.alphabet,
      node.is_root.then_some(&partition.gtr.pi),
    )?
  };

  if node.is_root {
    partition.nodes.get_mut(&node.key).unwrap().profile = msg_to_parent;
  } else {
    // Store msg_to_parent and propagate it through P(t)^T for the forward pass.
    let edge_key = get_exactly_one(&node.parent_edge_keys).expect("Only nodes with exactly one parent are supported");
    let branch_length = node.parent_edges[0].branch_length().unwrap_or(0.0);
    let branch_length = fix_branch_length(length, branch_length);
    let mut edge_data = partition.edges.remove(edge_key).unwrap();
    edge_data.msg_from_child = if partition.gtr.has_site_rates() {
      propagate_raw_per_site(
        &partition.gtr,
        branch_length,
        true,
        &msg_to_parent,
        edge_data.transmission.as_deref(),
      )
    } else {
      propagate_raw(
        &partition.gtr.expQt(branch_length).t().to_owned(),
        &msg_to_parent,
        edge_data.transmission.as_deref(),
      )
    };
    edge_data.msg_to_parent = msg_to_parent;
    partition.edges.insert(*edge_key, edge_data);
  }
  Ok(())
}

pub fn process_node_forward<N, E>(
  partition: &mut PartitionMarginalSparse,
  graph: &Graph<N, E, ()>,
  node: &GraphNodeForward<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  if !node.is_root {
    let mut variable_pos = btreemap! {};
    let mut ref_states: Vec<BTreeMap<usize, AsciiChar>> = vec![];
    let mut msgs_to_combine: Vec<MarginalSparseSeqDistribution> = vec![];
    let mut removed_edges = vec![];
    for (parent_key, edge_key) in &node.parent_keys {
      let mut parent_state: BTreeMap<usize, AsciiChar> = btreemap! {};
      let mut child_state: BTreeMap<usize, AsciiChar> = btreemap! {};
      let edge_data = partition.edges.remove(edge_key).unwrap();
      removed_edges.push((*edge_key, edge_data.clone()));

      let node_data = &partition.nodes[&node.key];
      // process substitutions first. these have defined states for both parent and child.
      for m in edge_data.fitch_subs() {
        let current_state = m.qry();
        variable_pos.insert(m.pos(), current_state);
        parent_state.entry(m.pos()).or_insert_with(|| m.reff());
        child_state.entry(m.pos()).or_insert(current_state);
      }
      // variable in parent
      for (pos, p) in &edge_data.msg_to_child.variable {
        // check if the node as a defined sequence at this position. If yes, state is same as parent unless there is a substitution (which is covered above)
        if !range_contains(&node_data.seq.non_char, *pos) {
          variable_pos.entry(*pos).or_insert(p.state);
          parent_state.entry(*pos).or_insert(p.state);
        }
      }
      // variable in child
      for (pos, p) in &edge_data.msg_to_parent.variable {
        variable_pos.entry(*pos).or_insert(p.state);
        child_state.entry(*pos).or_insert(p.state);
      }

      let edge_payload = graph.get_edge(*edge_key).unwrap().read_arc().payload().read_arc();
      let branch_length = edge_payload.branch_length().unwrap_or(0.0);

      msgs_to_combine.push(if partition.gtr.has_site_rates() {
        propagate_raw_per_site(
          &partition.gtr,
          branch_length,
          false,
          &edge_data.msg_to_child,
          edge_data.transmission.as_deref(),
        )
      } else {
        propagate_raw(
          &partition.gtr.expQt(branch_length),
          &edge_data.msg_to_child,
          edge_data.transmission.as_deref(),
        )
      });
      msgs_to_combine.push(edge_data.msg_to_parent.clone());

      ref_states.push(parent_state);
      ref_states.push(child_state);
    }

    for (edge_key, edge_data) in removed_edges {
      partition.edges.insert(edge_key, edge_data);
    }

    let composition = partition.nodes[&node.key].seq.composition.clone();
    let profile = combine_messages(
      &composition,
      &msgs_to_combine,
      &variable_pos,
      &ref_states,
      &partition.alphabet,
      None,
    )?;

    partition.nodes.get_mut(&node.key).unwrap().profile = profile;

    // Reconstruct MAP sequence for this node so resolve_map_state can read it.
    // Internal node seq.sequence is cleared after Fitch to save memory; restore it
    // from parent sequence + edge mutations + variable site MAP states.
    reconstruct_map_sequence(partition, node);

    // Compute ML subs for parent edges now that child profile is finalized
    for (parent_key, edge_key) in &node.parent_keys {
      let ml_subs = compute_ml_subs_for_edge(partition, *parent_key, node.key, *edge_key)?;
      if let Some(edge) = partition.edges.get_mut(edge_key) {
        edge.set_ml_subs(ml_subs);
      }
    }
  }

  // precalculate messages to children that summarize info from their siblings and the parent
  for child_edge_key in &node.child_edge_keys {
    let child_edge_data = partition.edges.remove(child_edge_key).unwrap();
    let seq_info = &partition.nodes[&node.key];
    let mut seq_dis = MarginalSparseSeqDistribution {
      variable: btreemap! {},
      variable_indel: btreemap! {},
      fixed: btreemap! {},
      fixed_counts: seq_info.seq.composition.clone(),
      log_lh: seq_info.profile.log_lh - child_edge_data.msg_from_child.log_lh,
    };

    let child_dis = child_edge_data.msg_from_child.clone();
    let mut parent_states: BTreeMap<usize, AsciiChar> = btreemap! {};
    let mut child_states: BTreeMap<usize, AsciiChar> = btreemap! {};
    let child_key = graph.get_target_node_key(*child_edge_key)?;
    let child_non_char = &partition.nodes[&child_key].seq.non_char;
    // handle substitutions first, these have defined states for both parent and child
    for sub in child_edge_data.fitch_subs() {
      child_states.insert(sub.pos(), sub.qry());
      parent_states.insert(sub.pos(), sub.reff());
    }
    // if not yet set, variable in parent is same as in child (no substitution);
    // skip positions where the child has no sequence (gap or N)
    for (pos, p) in &seq_info.profile.variable {
      if range_contains(child_non_char, *pos) {
        continue;
      }
      child_states.entry(*pos).or_insert(p.state);
      parent_states.entry(*pos).or_insert(p.state);
    }
    for (pos, p) in &child_dis.variable {
      if range_contains(child_non_char, *pos) {
        continue;
      }
      child_states.entry(*pos).or_insert(p.state);
      parent_states.entry(*pos).or_insert(p.state);
    }

    let mut delta_ll = 0.0;
    for (pos, pstate) in parent_states {
      let divisor = if let Some(dis) = child_dis.variable.get(&pos) {
        &dis.dis
      } else {
        &child_dis.fixed[&child_states[&pos]]
      };
      let numerator = if let Some(dis) = seq_info.profile.variable.get(&pos) {
        &dis.dis
      } else {
        &seq_info.profile.fixed[&pstate]
      };
      // Guard against zero divisor: clamp to smallest positive normal f64.
      // Matches the stabilization pattern in discrete.rs:206.
      let safe_divisor = divisor.mapv(|v| v.max(f64::MIN_POSITIVE));
      let dis = numerator / &safe_divisor;
      let norm = dis.sum();
      let dis = if norm > 0.0 && norm.is_finite() {
        delta_ll += norm.ln();
        dis / norm
      } else {
        delta_ll += f64::NEG_INFINITY;
        Array1::from_elem(dis.len(), 1.0 / dis.len() as f64)
      };
      seq_dis.variable.insert(pos, VarPos { dis, state: pstate });
      seq_dis.fixed_counts.adjust_count(pstate, -1);
    }
    for (s, p) in &seq_info.profile.fixed {
      let safe_fixed = child_dis.fixed[s].mapv(|v| v.max(f64::MIN_POSITIVE));
      let dis = p / &safe_fixed;
      let norm = dis.sum();
      let dis = if norm > 0.0 && norm.is_finite() {
        delta_ll += norm.ln() * (seq_dis.fixed_counts.get(*s).unwrap() as f64);
        dis / norm
      } else {
        delta_ll += f64::NEG_INFINITY;
        Array1::from_elem(dis.len(), 1.0 / dis.len() as f64)
      };
      seq_dis.fixed.insert(*s, dis);
    }
    seq_dis.log_lh += delta_ll;
    let mut updated_edge_data = child_edge_data;
    updated_edge_data.msg_to_child = seq_dis;
    partition.edges.insert(*child_edge_key, updated_edge_data);
  }
  Ok(())
}
