use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::partition::marginal::dense::partition::assign_sequence;
use crate::partition::marginal::shared::data::DenseInputs;
use crate::partition::marginal::shared::normalize::{
  forward_log_lh_add_normalization, forward_log_lh_remove_child, normalize_from_log, normalize_inplace,
};
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalForward};
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseLeafObs, DenseNodeState, DenseSeqDistribution,
  DenseSeqInfo,
};
use crate::seq::indel::{InDel, compute_node_ranges, resolve_indels_backward, resolve_indels_forward};
use eyre::Report;
use itertools::Itertools;
use ndarray::Array2;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPass, GraphPassBackwardContext, GraphPassForwardContext, GraphPassNodeOutput};
use treetime_primitives::LogLh;
use treetime_utils::interval::range_union::range_union;
use treetime_utils::make_internal_error;

#[expect(
  clippy::zero_sized_map_values,
  reason = "the graph pass takes a per-node input map; this pass has no per-node input and reads leaf observations by key"
)]
pub(crate) fn indexed_backward(
  inputs: &DenseInputs,
  gtr: &GTR,
  kind: IndexedKind<'_>,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<MarginalBackward<DenseNodeState, DenseEdgeBackward>, Report> {
  let min_branch_length = inputs.min_branch_length;
  let pass = GraphPass::new(graph)?;
  let outputs = match kind {
    IndexedKind::Dense {
      alphabet,
      length,
      leaves,
    } => pass.map_backward(
      &BTreeMap::new(),
      branch_lengths,
      |_| Ok(()),
      |context| {
        indexed_node_backward(
          gtr,
          min_branch_length,
          &context,
          |key| dense_leaf_backward(alphabet, &leaves[&key]),
          |children| backward_internal_dense(children, length),
        )
      },
    )?,
    IndexedKind::Discrete { leaves } => pass.map_backward(
      &BTreeMap::new(),
      branch_lengths,
      |_| Ok(()),
      |context| {
        indexed_node_backward(
          gtr,
          min_branch_length,
          &context,
          |key| Ok(discrete_leaf_backward(&leaves[&key])),
          |_| DenseSeqInfo::default(),
        )
      },
    )?,
  };
  Ok(MarginalBackward {
    node_states: outputs.nodes,
    backward: outputs.edges,
  })
}

#[expect(
  clippy::expect_used,
  reason = "the graph pass publishes every child edge message before its parent and gives every non-root node its parent edge"
)]
fn indexed_node_backward(
  gtr: &GTR,
  min_branch_length: f64,
  context: &GraphPassBackwardContext<'_, &(), f64, DenseNodeState, DenseEdgeBackward>,
  leaf_backward: impl Fn(GraphNodeKey) -> Result<(DenseNodeState, DenseSeqDistribution), Report>,
  internal_seq: impl Fn(&[&DenseNodeState]) -> DenseSeqInfo,
) -> Result<GraphPassNodeOutput<DenseNodeState, DenseEdgeBackward>, Report> {
  let (mut node, msg_to_parent) = if context.is_leaf {
    leaf_backward(context.key)?
  } else {
    let children = context.children.iter().map(|child| child.node).collect_vec();
    let seq = internal_seq(&children);

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
    let profile = DenseSeqDistribution { dis, log_lh };
    (
      DenseNodeState {
        seq,
        profile: profile.clone(),
      },
      profile,
    )
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

fn dense_leaf_backward(
  alphabet: &Alphabet,
  obs: &DenseLeafObs,
) -> Result<(DenseNodeState, DenseSeqDistribution), Report> {
  let message = DenseSeqDistribution {
    dis: alphabet.seq2prof(obs.sequence())?,
    log_lh: LogLh::ZERO,
  };
  let node = DenseNodeState {
    seq: obs.seq_info(),
    profile: DenseSeqDistribution::default(),
  };
  Ok((node, message))
}

fn discrete_leaf_backward(obs: &Array2<f64>) -> (DenseNodeState, DenseSeqDistribution) {
  let profile = DenseSeqDistribution::new(obs.clone(), LogLh::ZERO);
  let node = DenseNodeState {
    seq: DenseSeqInfo::default(),
    profile: profile.clone(),
  };
  (node, profile)
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

pub(crate) fn indexed_forward(
  inputs: &DenseInputs,
  gtr: &GTR,
  kind: IndexedKind<'_>,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
  backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
) -> Result<MarginalForward<DenseNodeState, DenseEdgeForward, DenseEdgeEstimate>, Report> {
  let min_branch_length = inputs.min_branch_length;
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_forward_owned(
    node_states,
    backward,
    |key| make_internal_error!("Partition node {key} is missing before the marginal forward pass"),
    |context| indexed_node_forward(gtr, min_branch_length, kind, branch_lengths, context),
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
  Ok(MarginalForward {
    node_states: outputs.nodes,
    forward,
    estimates,
  })
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn indexed_node_forward(
  gtr: &GTR,
  min_branch_length: f64,
  kind: IndexedKind<'_>,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
  context: GraphPassForwardContext<'_, DenseNodeState, DenseEdgeBackward, DenseNodeState>,
) -> Result<GraphPassNodeOutput<DenseNodeState, DenseEdgeForwardOut>, Report> {
  let mut node = context.input;

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

  if let IndexedKind::Dense { alphabet, .. } = kind {
    let indels = forward_post_dense(context.is_root, context.is_leaf, context.parent, &mut node, alphabet)?;
    if let Some((_, _, out)) = edge_out.as_mut() {
      out.indels = indels;
    }
  }

  let parent_message = edge_out.map(|(_, _, out)| out);
  Ok(GraphPassNodeOutput { node, parent_message })
}

#[derive(Clone, Copy, Debug)]
pub(crate) enum IndexedKind<'a> {
  Dense {
    alphabet: &'a Alphabet,
    length: usize,
    leaves: &'a BTreeMap<GraphNodeKey, DenseLeafObs>,
  },
  Discrete {
    leaves: &'a BTreeMap<GraphNodeKey, Array2<f64>>,
  },
}

struct DenseEdgeForwardOut {
  msg_to_child: DenseSeqDistribution,
  indels: Vec<InDel>,
}

#[allow(
  clippy::expect_used,
  reason = "expect on a value an upstream invariant guarantees is present"
)]
fn forward_post_dense(
  is_root: bool,
  is_leaf: bool,
  parent: Option<&DenseNodeState>,
  node: &mut DenseNodeState,
  alphabet: &Alphabet,
) -> Result<Vec<InDel>, Report> {
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
