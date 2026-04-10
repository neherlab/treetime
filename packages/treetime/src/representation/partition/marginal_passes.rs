use crate::hacks::fix_branch_length::fix_branch_length;
use crate::representation::partition::marginal_helpers::{EPS, combine_messages, propagate_raw, propagate_raw_per_site};
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::ExactStateCache;
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, VarPos};
use eyre::Report;
use maplit::btreemap;
use ndarray::Array1;
use std::collections::BTreeMap;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, Named};
use treetime_primitives::AsciiChar;
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range::range_contains;

pub fn process_node_backward<N, E>(
  partition: &mut PartitionMarginalSparse,
  graph: &Graph<N, E, ()>,
  node: &GraphNodeBackward<N, E, ()>,
  cache: &mut ExactStateCache,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let length = partition.length;
  let parent = if node.is_root {
    None
  } else {
    let edge_key = *get_exactly_one(&node.parent_edge_keys).expect("Only nodes with exactly one parent are supported");
    Some((graph.get_source_node_key(edge_key)?, edge_key))
  };

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
      .map(|(pos, p)| -> Result<_, Report> {
        let dis = alphabet.construct_profile(p.chars()).unwrap();
        let state = partition.node_canonical_state_known_parent(graph, node.key, parent, *pos, cache)?;
        Ok((*pos, VarPos { dis, state }))
      })
      .collect::<Result<_, _>>()?;

    MarginalSparseSeqDistribution {
      fixed_counts: node_data.seq.composition.clone(),
      variable,
      variable_indel: btreemap! {},
      fixed,
      log_lh: 0.0,
    }
  } else {
    // for internal nodes, combine the messages from the children.
    //
    // Spec: variable candidate positions come from
    //   (a) substitutions on edges to each child (node state = sub.reff),
    //   (b) positions variable in children (node state = child state when no sub on edge).
    //
    // `sub.reff()` and `msg_from_child.variable[pos].state` both carry the
    // current node's canonical state by construction, so there is no need to
    // walk the graph to derive it. `propagate_raw` preserves the `state` field
    // through edge transitions, so the child's state field in `msg_from_child`
    // equals the child's canonical state, which equals the node's state at a
    // position where the edge does not mutate.
    let mut variable_pos = btreemap! {};
    let mut child_states: Vec<BTreeMap<usize, AsciiChar>> = vec![];
    let mut child_messages: Vec<MarginalSparseSeqDistribution> = vec![];
    for (ci, (_child_key, edge_key)) in node.child_keys.iter().enumerate() {
      child_states.push(btreemap! {});
      let edge_data = &partition.edges[edge_key];
      for m in &edge_data.subs {
        variable_pos.insert(m.pos(), m.reff());
        child_states[ci].insert(m.pos(), m.qry());
      }
      for (pos, p) in &edge_data.msg_from_child.variable {
        variable_pos.entry(*pos).or_insert(p.state);
      }
      child_messages.push(edge_data.msg_from_child.clone());
    }

    // Resolve each child's state at each variable position.
    //
    // Precedence matches `child_canonical_state_from_parent_state()`:
    //   1. Child's gap range (highest)
    //   2. Child's unknown range
    //   3. Edge substitution (sub.qry)
    //   4. Parent (node) state
    //
    // The mask ranges take precedence over substitutions because a deletion or
    // unknown at the child masks the substitution model entirely: there is no
    // nucleotide at that position to receive the substituted state. Subs pre-
    // filled from the backward collection loop above are overwritten here when
    // a mask applies.
    let alphabet = &partition.alphabet;
    for (ci, (child_key, _)) in node.child_keys.iter().enumerate() {
      let states = &mut child_states[ci];
      let child_data = &partition.nodes[child_key];
      for (pos, parent_state) in &variable_pos {
        if range_contains(&child_data.seq.gaps, *pos) {
          states.insert(*pos, alphabet.gap());
        } else if range_contains(&child_data.seq.unknown, *pos) {
          states.insert(*pos, alphabet.unknown());
        } else if !states.contains_key(pos) {
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
    // data from children * gtr.pi as calculated above is the root profile.
    partition.nodes.get_mut(&node.key).unwrap().profile = msg_to_parent;
  } else {
    // what was calculated above is what is sent to the parent. we also calculate the propagated message to the parent (we need it in the forward pass).
    let edge_key = get_exactly_one(&node.parent_edge_keys).expect("Only nodes with exactly one parent are supported");
    let branch_length = node.parent_edges[0].branch_length().unwrap_or(0.0);
    let branch_length = fix_branch_length(length, branch_length);
    let edge_data = partition.edges.get_mut(edge_key).unwrap();
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
  }
  Ok(())
}

pub fn process_node_forward<N, E>(
  partition: &mut PartitionMarginalSparse,
  graph: &Graph<N, E, ()>,
  node: &GraphNodeForward<N, E, ()>,
  _cache: &mut ExactStateCache,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  if !node.is_root {
    // Combine parent-side propagated message with node's backward message.
    //
    // The `state` field on each `VarPos` is preserved by `propagate_raw` and
    // carries the node-side Fitch state by construction:
    //   * `msg_to_child.variable[pos].state` = parent's state at pos
    //   * `msg_to_parent.variable[pos].state` = node's state at pos
    //   * `edge.subs` defines both sides explicitly
    //
    // Spec: at positions without a sub on the edge, parent state == node state.
    // That makes the per-message reference states derivable from the message's
    // own state field plus the edge's sub list, without any graph traversal.
    let mut variable_pos = btreemap! {};
    let mut ref_states: Vec<BTreeMap<usize, AsciiChar>> = vec![];
    let mut msgs_to_combine: Vec<MarginalSparseSeqDistribution> = vec![];
    for (_parent_key, edge_key) in &node.parent_keys {
      let edge_data = &partition.edges[edge_key];
      let node_data = &partition.nodes[&node.key];

      let mut parent_state: BTreeMap<usize, AsciiChar> = btreemap! {};
      let mut child_state: BTreeMap<usize, AsciiChar> = btreemap! {};

      // Subs have both sides defined. Process first so entries win for later
      // sources that would otherwise only set one side.
      for m in &edge_data.subs {
        variable_pos.insert(m.pos(), m.qry());
        parent_state.entry(m.pos()).or_insert_with(|| m.reff());
        child_state.entry(m.pos()).or_insert_with(|| m.qry());
      }

      // Parent-side variable (msg_to_child). Skip positions where the node has
      // a gap or unknown: fixed-site lookups would fail at non-canonical chars.
      for (pos, p) in &edge_data.msg_to_child.variable {
        if range_contains(&node_data.seq.non_char, *pos) {
          continue;
        }
        variable_pos.entry(*pos).or_insert(p.state);
        parent_state.entry(*pos).or_insert(p.state);
      }

      // Node-side variable (msg_to_parent): node state = p.state.
      for (pos, p) in &edge_data.msg_to_parent.variable {
        variable_pos.entry(*pos).or_insert(p.state);
        child_state.entry(*pos).or_insert(p.state);
      }

      // Fill cross-side fallbacks: at positions with no sub on this edge, the
      // parent state and child state coincide with the message's own state.
      for (pos, state) in &variable_pos {
        parent_state.entry(*pos).or_insert(*state);
        child_state.entry(*pos).or_insert(*state);
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
  }

  // precalculate messages to children that summarize info from their siblings and the parent.
  //
  // `parent_states` = current node's state at each variable position (this node is
  // the parent of the outgoing edge). `child_states` = child's state at the same
  // positions. Both come from the existing message `state` fields and the edge
  // `subs` list; no graph traversal is needed.
  for child_edge_key in &node.child_edge_keys {
    let child_edge_data = partition.edges[child_edge_key].clone();
    let seq_info = &partition.nodes[&node.key];
    let mut seq_dis = MarginalSparseSeqDistribution {
      variable: btreemap! {},
      variable_indel: btreemap! {},
      fixed: btreemap! {},
      fixed_counts: seq_info.seq.composition.clone(),
      log_lh: seq_info.profile.log_lh - child_edge_data.msg_from_child.log_lh,
    };

    let child_dis = child_edge_data.msg_from_child.clone();
    let child_data = &partition.nodes[&graph.get_target_node_key(*child_edge_key)?];
    let alphabet = &partition.alphabet;

    let mut parent_states: BTreeMap<usize, AsciiChar> = btreemap! {};
    let mut child_states: BTreeMap<usize, AsciiChar> = btreemap! {};

    // Subs define both sides.
    for sub in &child_edge_data.subs {
      parent_states.insert(sub.pos(), sub.reff());
      child_states.insert(sub.pos(), sub.qry());
    }
    // Current node's variable: parent state = p.state; child state equals it
    // unless overridden by a sub or by the child's gap/unknown mask.
    for (pos, p) in &seq_info.profile.variable {
      parent_states.entry(*pos).or_insert(p.state);
      child_states.entry(*pos).or_insert(p.state);
    }
    // Child's upward variable: child state = p.state; parent state equals it
    // when the edge does not mutate at this position.
    for (pos, p) in &child_dis.variable {
      child_states.entry(*pos).or_insert(p.state);
      parent_states.entry(*pos).or_insert(p.state);
    }
    // Apply child gap/unknown masks unconditionally. Precedence matches
    // `child_canonical_state_from_parent_state()`: gap > unknown > sub > p.state.
    // Masks override substitutions because a deletion or unknown at the child
    // masks the substitution model entirely - the nucleotide does not exist
    // there to receive the substituted state.
    for (pos, state) in &mut child_states {
      if range_contains(&child_data.seq.gaps, *pos) {
        *state = alphabet.gap();
      } else if range_contains(&child_data.seq.unknown, *pos) {
        *state = alphabet.unknown();
      }
    }

    let neutral = Array1::from_elem(partition.alphabet.n_canonical(), 1.0);

    let mut delta_ll = 0.0;
    for (pos, pstate) in parent_states {
      let divisor = if let Some(dis) = child_dis.variable.get(&pos) {
        &dis.dis
      } else if partition.alphabet.is_canonical(child_states[&pos]) {
        &child_dis.fixed[&child_states[&pos]]
      } else {
        &neutral
      };
      let numerator = if let Some(dis) = seq_info.profile.variable.get(&pos) {
        &dis.dis
      } else if partition.alphabet.is_canonical(pstate) {
        &seq_info.profile.fixed[&pstate]
      } else {
        &neutral
      };
      let dis = numerator / divisor;
      let norm = dis.sum();
      delta_ll += norm.ln();
      let dis = dis / norm;
      let max_prob = dis.iter().copied().fold(0.0_f64, f64::max);
      let map_state = partition
        .alphabet
        .char(treetime_utils::array::ndarray::argmax_first(&dis.view()).unwrap_or(0));
      if max_prob < (1.0 - EPS) || map_state != pstate {
        seq_dis.variable.insert(pos, VarPos { dis, state: pstate });
      }
      seq_dis.fixed_counts.adjust_count(pstate, -1);
    }
    for (s, p) in &seq_info.profile.fixed {
      let dis = p / &child_dis.fixed[s];
      let norm = dis.sum();
      delta_ll += norm.ln() * (seq_dis.fixed_counts.get(*s).unwrap() as f64);
      seq_dis.fixed.insert(*s, dis / norm);
    }
    seq_dis.log_lh += delta_ll;
    partition.edges.get_mut(child_edge_key).unwrap().msg_to_child = seq_dis;
  }
  Ok(())
}
