use crate::hacks::fix_branch_length::fix_branch_length;
use crate::representation::partition::marginal_helpers::{combine_messages, propagate_raw, propagate_raw_per_site};
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::ExactStateCache;
use crate::representation::payload::sparse::{MarginalSparseSeqDistribution, VarPos};
use eyre::Report;
use maplit::btreemap;
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
  _graph: &Graph<N, E, ()>,
  node: &GraphNodeBackward<N, E, ()>,
  _cache: &ExactStateCache,
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
        let state = p.get_one();
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
    // for internal nodes, combine the messages from the children
    let mut variable_pos = btreemap! {};
    let mut child_states = vec![];
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

    let alphabet = &partition.alphabet;
    for (ci, (child_key, _)) in node.child_keys.iter().enumerate() {
      let states = &mut child_states[ci];
      let child_data = &partition.nodes[child_key];
      for pos in variable_pos.keys() {
        if !states.contains_key(pos) && range_contains(&child_data.seq.non_char, *pos) {
          if range_contains(&child_data.seq.gaps, *pos) {
            states.insert(*pos, alphabet.gap());
          } else {
            states.insert(*pos, alphabet.unknown());
          }
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
  _cache: &ExactStateCache,
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
    for (_, edge_key) in &node.parent_keys {
      let mut parent_state: BTreeMap<usize, AsciiChar> = btreemap! {};
      let mut child_state: BTreeMap<usize, AsciiChar> = btreemap! {};
      let edge_data = partition.edges.remove(edge_key).unwrap();
      removed_edges.push((*edge_key, edge_data.clone()));

      let node_data = &partition.nodes[&node.key];
      for (pos, p) in &edge_data.msg_to_child.variable {
        if !range_contains(&node_data.seq.gaps, *pos) {
          variable_pos.entry(*pos).or_insert(p.state);
          parent_state.insert(*pos, p.state);
        }
      }
      for m in &edge_data.subs {
        variable_pos.insert(m.pos(), m.qry());
        parent_state.entry(m.pos()).or_insert_with(|| m.reff());
        child_state.entry(m.pos()).or_insert_with(|| m.qry());
      }
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
    for (pos, p) in &seq_info.profile.variable {
      child_states.insert(*pos, p.state);
      parent_states.insert(*pos, p.state);
    }
    for (pos, p) in &child_dis.variable {
      child_states.insert(*pos, p.state);
      parent_states.entry(*pos).or_insert(p.state);
    }
    for sub in &child_edge_data.subs {
      child_states.insert(sub.pos(), sub.qry());
      parent_states.insert(sub.pos(), sub.reff());
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
      let dis = numerator / divisor;
      let norm = dis.sum();
      delta_ll += norm.ln();
      seq_dis.variable.insert(
        pos,
        VarPos {
          dis: dis / norm,
          state: pstate,
        },
      );
      seq_dis.fixed_counts.adjust_count(pstate, -1);
    }
    for (s, p) in &seq_info.profile.fixed {
      let dis = p / &child_dis.fixed[s];
      let norm = dis.sum();
      delta_ll += norm.ln() * (seq_dis.fixed_counts.get(*s).unwrap() as f64);
      seq_dis.fixed.insert(*s, dis / norm);
    }
    seq_dis.log_lh += delta_ll;
    let mut updated_edge_data = child_edge_data;
    updated_edge_data.msg_to_child = seq_dis;
    partition.edges.insert(*child_edge_key, updated_edge_data);
  }
  Ok(())
}
