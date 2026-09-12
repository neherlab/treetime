use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch_indel::{compute_node_ranges, resolve_indels_backward, resolve_indels_forward};
use crate::gtr::gtr::GTR;
use crate::partition::marginal::dense::partition::assign_sequence;
use crate::partition::marginal::shared::data::DenseInputs;
use crate::partition::marginal::shared::normalize::{
  forward_log_lh_add_normalization, forward_log_lh_remove_child, normalize_from_log, normalize_inplace,
};
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseNodeState, DenseSeqDistribution, DenseSeqInfo,
};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPass, GraphPassBackwardContext, GraphPassForwardContext, GraphPassNodeOutput};
use treetime_primitives::LogLh;
use treetime_utils::interval::range_union::range_union;

/// The two representations that ride the shared indexed (full-profile) marginal machinery. Dense
/// carries sequences, gaps, and indels; discrete carries a single categorical site and no indel work.
/// Representation is dispatched once, at the driver entry, and branches only at the three points where
/// the two genuinely differ (leaf profile source, backward site-info, forward post-processing).
#[derive(Clone, Copy, Debug)]
pub enum IndexedKind {
  Dense,
  Discrete,
}

/// Run the indexed marginal backward pass over borrowed inputs and node states, returning the updated
/// node states and the per-edge backward messages as distinct owned values.
pub fn indexed_backward(
  inputs: &DenseInputs,
  alphabet: &Alphabet,
  length: usize,
  kind: IndexedKind,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
) -> Result<(BTreeMap<GraphNodeKey, DenseNodeState>, BTreeMap<GraphEdgeKey, DenseEdgeBackward>), Report> {
  let gtr = &inputs.gtr;
  let min_branch_length = inputs.min_branch_length;
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_backward(
    node_states,
    branch_lengths,
    |_| Ok(DenseNodeState::empty()),
    |context| indexed_node_backward(gtr, min_branch_length, alphabet, length, kind, &context),
  )?;
  Ok((outputs.nodes, outputs.edges))
}

fn indexed_node_backward(
  gtr: &GTR,
  min_branch_length: f64,
  alphabet: &Alphabet,
  length: usize,
  kind: IndexedKind,
  context: &GraphPassBackwardContext<'_, DenseNodeState, f64, DenseNodeState, DenseEdgeBackward>,
) -> Result<GraphPassNodeOutput<DenseNodeState, DenseEdgeBackward>, Report> {
  let mut node = context.input.clone();
  let msg_to_parent = if context.is_leaf {
    match kind {
      IndexedKind::Dense => DenseSeqDistribution {
        dis: alphabet.seq2prof(&node.seq.sequence)?,
        log_lh: LogLh::ZERO,
      },
      IndexedKind::Discrete => node.profile.clone(),
    }
  } else {
    node.seq = match kind {
      IndexedKind::Dense => {
        let children = context.children.iter().map(|child| child.node).collect_vec();
        backward_internal_dense(&children, length)
      },
      IndexedKind::Discrete => DenseSeqInfo::default(),
    };

    // Children arrive in the graph's canonical `children_of` order, so the per-child log-space product
    // folds in the same order regardless of thread count, keeping the result byte-for-byte identical.
    let child_edges = context
      .children
      .iter()
      .map(|child| {
        child
          .edge
          .expect("Backward child edge message must be published before its parent")
      })
      .collect_vec();
    let first_edge = child_edges.first().expect("Internal node must have children");
    let mut log_dis = first_edge.msg_from_child.dis.mapv(f64::ln);
    for edge in child_edges.iter().skip(1) {
      log_dis += &edge.msg_from_child.dis.mapv(f64::ln);
    }
    let (dis, delta_ll) = normalize_from_log(&log_dis);
    let log_lh = child_edges.iter().map(|edge| edge.msg_from_child.log_lh).sum::<LogLh>() + LogLh::new(delta_ll);
    node.profile = DenseSeqDistribution {
      dis: dis.clone(),
      log_lh,
    };
    DenseSeqDistribution { dis, log_lh }
  };

  let parent_message = if context.is_root {
    let mut dis = &msg_to_parent.dis * &gtr.pi;
    let delta_ll = normalize_inplace(&mut dis);
    node.profile = DenseSeqDistribution {
      dis,
      log_lh: msg_to_parent.log_lh + LogLh::new(delta_ll),
    };
    None
  } else {
    let (_edge_key, branch_length) = context.parent_edge.expect("Non-root node must own its parent edge");
    let branch_length = branch_length.max(min_branch_length);
    Some(DenseEdgeBackward {
      msg_from_child: DenseSeqDistribution {
        dis: gtr.propagate_profile(&msg_to_parent.dis, branch_length, false),
        log_lh: msg_to_parent.log_lh,
      },
      msg_to_parent,
    })
  };

  Ok(GraphPassNodeOutput { node, parent_message })
}

fn backward_internal_dense(children: &[&DenseNodeState], length: usize) -> DenseSeqInfo {
  let child_non_chars = children.iter().map(|child| &child.seq.non_char).collect_vec();
  let child_gaps = children.iter().map(|child| &child.seq.gaps).collect_vec();
  let ranges = compute_node_ranges(&child_non_chars, &child_gaps);
  let child_unknown = children.iter().map(|child| &child.seq.unknown).collect_vec();
  let child_variable_indels = children.iter().map(|child| &child.seq.variable_indel).collect_vec();
  let indels = resolve_indels_backward(&child_gaps, &child_unknown, &child_variable_indels, length);
  DenseSeqInfo {
    gaps: indels.resolved_gaps,
    unknown: ranges.unknown,
    non_char: ranges.non_char,
    variable_indel: indels.variable_indel,
    ..DenseSeqInfo::default()
  }
}

/// Run the indexed marginal forward pass over borrowed inputs, node states, and backward messages,
/// returning the updated node states, the per-edge forward messages, and the per-edge estimates
/// (indels) as distinct owned values.
pub fn indexed_forward(
  inputs: &DenseInputs,
  alphabet: &Alphabet,
  kind: IndexedKind,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
  backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
) -> Result<
  (
    BTreeMap<GraphNodeKey, DenseNodeState>,
    BTreeMap<GraphEdgeKey, DenseEdgeForward>,
    BTreeMap<GraphEdgeKey, DenseEdgeEstimate>,
  ),
  Report,
> {
  let gtr = &inputs.gtr;
  let min_branch_length = inputs.min_branch_length;
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_forward(
    node_states,
    backward,
    |key| treetime_utils::make_internal_error!("Partition node {key} is missing before the marginal forward pass"),
    |context| indexed_node_forward(gtr, min_branch_length, alphabet, kind, branch_lengths, &context),
  )?;

  let mut forward = BTreeMap::new();
  let mut estimates = BTreeMap::new();
  for (edge_key, out) in outputs.edges {
    forward.insert(
      edge_key,
      DenseEdgeForward {
        msg_to_child: out.msg_to_child,
      },
    );
    estimates.insert(edge_key, DenseEdgeEstimate { indels: out.indels });
  }
  Ok((outputs.nodes, forward, estimates))
}

/// Combined per-edge output of the forward node visit, split by the driver into the distinct
/// forward-message and estimate owners.
struct DenseEdgeForwardOut {
  msg_to_child: DenseSeqDistribution,
  indels: Vec<crate::seq::indel::InDel>,
}

fn indexed_node_forward(
  gtr: &GTR,
  min_branch_length: f64,
  alphabet: &Alphabet,
  kind: IndexedKind,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: &GraphPassForwardContext<'_, DenseNodeState, DenseEdgeBackward, DenseNodeState>,
) -> Result<GraphPassNodeOutput<DenseNodeState, DenseEdgeForwardOut>, Report> {
  let mut node = context.input.clone();

  let mut edge_out = context.parent_edge.map(|(edge_key, backward)| {
    (
      edge_key,
      backward,
      DenseEdgeForwardOut {
        msg_to_child: DenseSeqDistribution::default(),
        indels: vec![],
      },
    )
  });

  if let Some((edge_key, backward, out)) = edge_out.as_mut() {
    let parent = context.parent.expect("Non-root node must have a parent");
    let safe_child = backward.msg_from_child.dis.mapv(|value| value.max(f64::MIN_POSITIVE));
    let mut dis = &parent.profile.dis / &safe_child;
    let delta_ll = normalize_inplace(&mut dis);
    let log_lh = forward_log_lh_remove_child(parent.profile.log_lh, backward.msg_from_child.log_lh);
    let log_lh = forward_log_lh_add_normalization(log_lh, delta_ll);
    out.msg_to_child = DenseSeqDistribution { dis, log_lh };
    let branch_length = branch_lengths[&*edge_key].max(min_branch_length);
    let msg_child = gtr.evolve(&out.msg_to_child.dis, branch_length, false);
    let mut dis = &backward.msg_to_parent.dis * &msg_child;
    let delta_ll = normalize_inplace(&mut dis);
    node.profile = DenseSeqDistribution {
      dis,
      log_lh: backward.msg_to_parent.log_lh + out.msg_to_child.log_lh + LogLh::new(delta_ll),
    };
  }

  if let IndexedKind::Dense = kind {
    let indels = forward_post_dense(context.is_root, context.is_leaf, context.parent, &mut node, alphabet)?;
    if let Some((_, _, out)) = edge_out.as_mut() {
      out.indels = indels;
    }
  }

  let parent_message = edge_out.map(|(_, _, out)| out);
  Ok(GraphPassNodeOutput { node, parent_message })
}

/// Dense forward post-processing: reconstruct the node sequence from the parent's reconstructed
/// sequence and gap structure and derive the branch indels. This keeps the per-node parent-sequence and
/// gap dependency intact; indel work runs during the forward node visit, not after the whole pass.
fn forward_post_dense(
  is_root: bool,
  is_leaf: bool,
  parent: Option<&DenseNodeState>,
  node: &mut DenseNodeState,
  alphabet: &Alphabet,
) -> Result<Vec<crate::seq::indel::InDel>, Report> {
  if is_root {
    node.seq.variable_indel.clear();
    node.seq.sequence = assign_sequence(node, alphabet);
    node.seq.non_char = range_union(&[node.seq.gaps.clone(), node.seq.unknown.clone()]);
    return Ok(vec![]);
  }

  if !is_leaf {
    node.seq.sequence = assign_sequence(node, alphabet);
  }
  let parent = parent.expect("Non-root dense node must have a parent");
  let variable_indel = std::mem::take(&mut node.seq.variable_indel);
  let (indels, new_gaps) = resolve_indels_forward(
    &variable_indel,
    &node.seq.gaps,
    &node.seq.non_char,
    &parent.seq.gaps,
    &parent.seq.sequence,
    &node.seq.sequence,
  );
  node.seq.gaps = new_gaps;
  node.seq.non_char = range_union(&[node.seq.gaps.clone(), node.seq.unknown.clone()]);
  for gap in &node.seq.gaps {
    node.seq.sequence[gap.0..gap.1].fill(alphabet.gap());
  }
  for unknown in &node.seq.unknown {
    node.seq.sequence[unknown.0..unknown.1].fill(alphabet.unknown());
  }
  Ok(indels)
}
