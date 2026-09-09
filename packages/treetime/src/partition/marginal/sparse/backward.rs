use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::partition::marginal::sparse::message::{combine_messages, propagate_raw, propagate_raw_per_site};
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::sparse::{SparseEdgePartition, SparseNodePartition, SparseSeqDistribution, VarPos};
use eyre::Report;
use maplit::btreemap;
use std::collections::{BTreeMap, BTreeSet};
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_graph::pass::{GraphPass, GraphPassBackwardContext, GraphPassNodeOutput};
use treetime_primitives::LogLh;
use treetime_utils::interval::range::range_contains;

pub fn process_backward_indexed<N, E>(
  partition: &mut PartitionMarginalSparse,
  graph: &Graph<N, E, ()>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let alphabet = partition.alphabet.clone();
  let gtr = partition.gtr.clone();
  let length = partition.length;
  let (nodes, edges) = (&mut partition.nodes, &mut partition.edges);
  let pass = GraphPass::new(graph, nodes, edges, |key| {
    treetime_utils::make_internal_error!("Partition node {key} is missing before the sparse marginal pass")
  })?;
  let outputs = pass.try_map_backward(|context| {
    process_node_backward_indexed(graph, &alphabet, &gtr, length, branch_lengths, context)
  })?;
  partition.nodes = outputs.nodes;
  partition.edges = outputs.edges;
  Ok(())
}

fn process_node_backward_indexed<N, E>(
  graph: &Graph<N, E, ()>,
  alphabet: &Alphabet,
  gtr: &GTR,
  length: usize,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: GraphPassBackwardContext<
    '_,
    SparseNodePartition,
    SparseEdgePartition,
    SparseNodePartition,
    SparseEdgePartition,
  >,
) -> Result<GraphPassNodeOutput<SparseNodePartition, SparseEdgePartition>, Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let mut node = context.input;
  let graph_node = graph.get_node(context.key).expect("Indexed node must exist in graph");
  let graph_node = graph_node.read_arc();
  let msg_to_parent = if graph_node.is_leaf() {
    let fixed = alphabet
      .determined()
      .map(|state| Ok((state, alphabet.get_profile(state)?.clone())))
      .collect::<Result<_, Report>>()?;
    let variable = node
      .seq
      .fitch
      .variable
      .iter()
      .map(|(pos, profile)| {
        let dis = alphabet.construct_profile(profile.chars()).unwrap();
        let state = node
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
      fixed_counts: node.seq.composition.clone(),
      variable,
      variable_indel: BTreeSet::new(),
      fixed,
      log_lh: LogLh::ZERO,
    }
  } else {
    let mut variable_pos = btreemap! {};
    let child_pairs = graph.children_of(&graph_node);

    // The value engine hands the completed children in its own topology order, which may differ from
    // `children_of`. Index them by key so every child is fetched, and folded, in the same canonical
    // `children_of` order as before, keeping the floating-point result byte-for-byte identical.
    let child_nodes: BTreeMap<_, _> = context
      .children
      .iter()
      .map(|child| (child.node_key, child.node))
      .collect();
    let child_edge_messages: BTreeMap<_, _> = context
      .children
      .iter()
      .filter_map(|child| child.edge.map(|edge| (child.edge_key, edge)))
      .collect();

    let mut child_states = vec![btreemap! {}; child_pairs.len()];
    let mut child_messages = Vec::with_capacity(child_pairs.len());
    let mut child_keys = Vec::with_capacity(child_pairs.len());

    for (ci, (child, edge)) in child_pairs.iter().enumerate() {
      let child_key = child.read_arc().key();
      let edge_key = edge.read_arc().key();
      let edge_data = *child_edge_messages
        .get(&edge_key)
        .expect("Backward child edge message must be published before its parent");
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
      let child_data = *child_nodes
        .get(child_key)
        .expect("Backward child node output must be published before its parent");
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
      &node.seq.composition,
      &child_messages,
      &variable_pos,
      &child_states,
      alphabet,
      graph_node.is_root().then_some(&gtr.pi),
    )?
  };

  let parent_message = if graph_node.is_root() {
    node.profile = msg_to_parent;
    None
  } else {
    // Reuse this node's moved-in parent-edge input and overwrite only the two message fields, exactly
    // as the in-place engine did. The edge already holds `fitch_subs` and `transmission` written by the
    // Fitch pre-pass, plus other fields, and a fresh edge would destroy them and corrupt the result.
    let (edge_key, mut edge_data) = context.parent_edge.expect("Non-root node must own its parent edge");
    let branch_length = fix_branch_length(length, branch_lengths[&edge_key]);
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
    Some(edge_data)
  };

  Ok(GraphPassNodeOutput { node, parent_message })
}
