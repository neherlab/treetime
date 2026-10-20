use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::partition::marginal::dense::partition::assign_sequence;
use crate::partition::marginal::shared::data::{DenseInputs, count_transitions_dense};
use crate::partition::marginal::shared::normalize::{
  divide_out, forward_log_lh_add_normalization, forward_log_lh_remove_child, normalize_from_log, normalize_inplace,
};
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalForward, MarginalPasses};
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseLeafObs, DenseNodeState, DenseSeqDistribution,
  DenseSeqInfo,
};
use crate::seq::indel::{InDel, compute_node_ranges, resolve_indels_backward, resolve_indels_forward};
use eyre::Report;
use itertools::{Itertools, izip};
use ndarray::Array2;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::NodeEdgeMaps;
use treetime_primitives::LogLh;
use treetime_utils::interval::range_union::range_union;
use treetime_utils::{make_internal_error, make_internal_report};

impl<P: IndexedPartition> MarginalPasses for P {
  type Node = DenseNodeState;
  type BackwardInput = ();
  type Backward = DenseEdgeBackward;
  type Forward = DenseEdgeForward;
  type Estimate = DenseEdgeEstimate;

  fn backward_input(_node_states: &BTreeMap<GraphNodeKey, DenseNodeState>) -> &Self::BackwardInput {
    &()
  }

  fn fresh_backward_input(_node_states: &BTreeMap<GraphNodeKey, DenseNodeState>) -> Self::BackwardInput {}

  fn marginal_backward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    (): &(),
  ) -> Result<MarginalBackward<DenseNodeState, DenseEdgeBackward>, Report> {
    let kind = self.indexed_kind();
    gtr.check_site_count(kind.site_count())?;
    indexed_backward(self.inputs(), gtr, kind, graph, branch_lengths)
  }

  fn marginal_forward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
  ) -> Result<MarginalForward<DenseNodeState, DenseEdgeForward, DenseEdgeEstimate>, Report> {
    let kind = self.indexed_kind();
    gtr.check_site_count(kind.site_count())?;
    indexed_forward(self.inputs(), gtr, kind, graph, branch_lengths, node_states, backward)
  }

  fn count_transitions(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
    node_states: &BTreeMap<GraphNodeKey, DenseNodeState>,
    backward: &BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
    forward: &BTreeMap<GraphEdgeKey, DenseEdgeForward>,
  ) -> Result<MutationCounts, Report> {
    count_transitions_dense(
      self.inputs(),
      gtr,
      graph,
      branch_lengths,
      node_states,
      backward,
      forward,
    )
  }
}

pub(crate) fn indexed_backward(
  inputs: &DenseInputs,
  gtr: &GTR,
  kind: IndexedKind<'_>,
  graph: &Graph,
  branch_lengths: &BTreeMap<GraphEdgeKey, f64>,
) -> Result<MarginalBackward<DenseNodeState, DenseEdgeBackward>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let pass = IndexedBackwardPass {
    schedule: &schedule,
    gtr,
    min_branch_length: inputs.min_branch_length,
    branch_lengths,
  };
  let outputs = match kind {
    IndexedKind::Dense {
      alphabet,
      length,
      leaves,
    } => schedule.scan_up(leaf_inputs(&schedule, leaves), |key, leaf, children| {
      pass.node(
        key,
        leaf,
        &children,
        |obs| dense_leaf_backward(alphabet, obs),
        |children| backward_internal_dense(children, length),
      )
    })?,
    IndexedKind::Discrete { leaves } => schedule.scan_up(leaf_inputs(&schedule, leaves), |key, leaf, children| {
      pass.node(
        key,
        leaf,
        &children,
        |obs| Ok(discrete_leaf_backward(obs)),
        |_| DenseSeqInfo::default(),
      )
    })?,
  };
  let NodeEdgeMaps { nodes, edges } = outputs.unwrap_shared().split_parent_edges()?;
  Ok(MarginalBackward {
    node_states: nodes,
    backward: edges,
  })
}

fn leaf_inputs<'a, L>(
  schedule: &'a GraphSchedule<'_>,
  leaves: &'a BTreeMap<GraphNodeKey, L>,
) -> impl Iterator<Item = (GraphNodeKey, Option<&'a L>)> {
  schedule.preorder().iter().map(|&key| (key, leaves.get(&key)))
}

type IndexedBackwardOutput = Arc<(DenseNodeState, Vec<DenseEdgeBackward>)>;

struct IndexedBackwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  gtr: &'a GTR,
  min_branch_length: f64,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, f64>,
}

impl IndexedBackwardPass<'_> {
  fn node<L>(
    &self,
    key: GraphNodeKey,
    leaf: Option<&L>,
    children: &[IndexedBackwardOutput],
    leaf_backward: impl Fn(&L) -> Result<(DenseNodeState, DenseSeqDistribution), Report>,
    internal_seq: impl Fn(&[&DenseNodeState]) -> DenseSeqInfo,
  ) -> Result<(IndexedBackwardOutput, IndexedBackwardOutput), Report> {
    let (mut node, msg_to_parent) = if self.schedule.children(key).is_empty() {
      let leaf =
        leaf.ok_or_else(|| make_internal_report!("Leaf {key} has no observation before the marginal backward pass"))?;
      leaf_backward(leaf)?
    } else {
      let child_nodes = children.iter().map(|child| &child.0).collect_vec();
      let seq = internal_seq(&child_nodes);
      let child_edges = izip!(self.schedule.children(key), children)
        .map(|(&(child_key, edge_key), child)| {
          child
            .1
            .get(self.schedule.parent_slot(edge_key))
            .ok_or_else(|| make_internal_report!("Backward child {child_key} published no message for its parent edge"))
        })
        .collect::<Result<Vec<_>, Report>>()?;
      let Some((first_edge, other_edges)) = child_edges.split_first() else {
        return make_internal_error!("Internal node {key} has no children");
      };
      let mut log_dis = first_edge.msg_from_child.dis.mapv(f64::ln);
      for edge in other_edges {
        log_dis.zip_mut_with(&edge.msg_from_child.dis, |l, &p| *l += p.ln());
      }
      let (dis, delta_ll) = normalize_from_log(log_dis);
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

    let parent_message = if let Some(&(_, edge_key)) = self.schedule.parents(key).first() {
      let branch_length = self.branch_lengths[&edge_key].max(self.min_branch_length);
      vec![DenseEdgeBackward {
        msg_from_child: DenseSeqDistribution {
          dis: self.gtr.propagate_profile(&msg_to_parent.dis, branch_length, false),
          log_lh: msg_to_parent.log_lh,
        },
        msg_to_parent,
      }]
    } else {
      let mut dis = &msg_to_parent.dis * self.gtr.pi();
      let delta_ll = normalize_inplace(&mut dis);
      node.profile = DenseSeqDistribution {
        dis,
        log_lh: msg_to_parent.log_lh + LogLh::new(delta_ll),
      };
      vec![]
    };

    let output = Arc::new((node, parent_message));
    Ok((Arc::clone(&output), output))
  }
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
  let schedule = GraphSchedule::new(graph)?;
  schedule.require_tree()?;
  let pass = IndexedForwardPass {
    schedule: &schedule,
    gtr,
    min_branch_length: inputs.min_branch_length,
    kind,
    branch_lengths,
  };
  let backward = backward.iter().map(|(&edge_key, edge)| (edge_key, edge)).collect();
  let outputs = schedule.scan_down(
    schedule.join_parent_edges(node_states, backward)?,
    |key, (node, mut backward), parents: Vec<IndexedForwardOutput>| {
      pass.node(key, node, backward.pop(), parents.first().map(|parent| &parent.0))
    },
  )?;
  let NodeEdgeMaps { nodes, edges } = outputs.unwrap_shared().split_parent_edges()?;

  let mut forward = BTreeMap::new();
  let mut estimates = BTreeMap::new();
  for (edge_key, out) in edges {
    forward.insert(
      edge_key,
      DenseEdgeForward {
        msg_to_child: out.msg_to_child,
      },
    );
    estimates.insert(edge_key, DenseEdgeEstimate { indels: out.indels });
  }
  Ok(MarginalForward {
    node_states: nodes,
    forward,
    estimates,
  })
}

type IndexedForwardOutput = Arc<(DenseNodeState, Vec<DenseEdgeForwardOut>)>;

struct IndexedForwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  gtr: &'a GTR,
  min_branch_length: f64,
  kind: IndexedKind<'a>,
  branch_lengths: &'a BTreeMap<GraphEdgeKey, f64>,
}

impl IndexedForwardPass<'_> {
  fn node(
    &self,
    key: GraphNodeKey,
    mut node: DenseNodeState,
    backward: Option<&DenseEdgeBackward>,
    parent: Option<&DenseNodeState>,
  ) -> Result<(IndexedForwardOutput, IndexedForwardOutput), Report> {
    let msg_to_child = match (self.schedule.parents(key).first(), backward, parent) {
      (Some(&(_, edge_key)), Some(backward), Some(parent)) => {
        let mut dis = divide_out(&parent.profile.dis, &backward.msg_from_child.dis);
        let delta_ll = normalize_inplace(&mut dis);
        let log_lh = forward_log_lh_remove_child(parent.profile.log_lh, backward.msg_from_child.log_lh);
        let log_lh = forward_log_lh_add_normalization(log_lh, delta_ll);
        let msg_to_child = DenseSeqDistribution { dis, log_lh };
        let branch_length = self.branch_lengths[&edge_key].max(self.min_branch_length);
        let msg_child = self.gtr.evolve(&msg_to_child.dis, branch_length, false);
        let mut dis = &backward.msg_to_parent.dis * &msg_child;
        let delta_ll = normalize_inplace(&mut dis);
        node.profile = DenseSeqDistribution {
          dis,
          log_lh: backward.msg_to_parent.log_lh + msg_to_child.log_lh + LogLh::new(delta_ll),
        };
        Some(msg_to_child)
      },
      (None, None, None) => None,
      _ => {
        return Err(make_internal_report!(
          "Marginal forward pass received an inconsistent parent edge and parent message at node {key}"
        ));
      },
    };

    let indels = match self.kind {
      IndexedKind::Dense { alphabet, .. } => {
        forward_post_dense(self.schedule.children(key).is_empty(), parent, &mut node, alphabet)?
      },
      IndexedKind::Discrete { .. } => vec![],
    };

    let edges = msg_to_child
      .map(|msg_to_child| DenseEdgeForwardOut { msg_to_child, indels })
      .into_iter()
      .collect_vec();
    let output = Arc::new((node, edges));
    Ok((Arc::clone(&output), output))
  }
}

struct DenseEdgeForwardOut {
  msg_to_child: DenseSeqDistribution,
  indels: Vec<InDel>,
}

fn forward_post_dense(
  is_leaf: bool,
  parent: Option<&DenseNodeState>,
  node: &mut DenseNodeState,
  alphabet: &Alphabet,
) -> Result<Vec<InDel>, Report> {
  let Some(parent) = parent else {
    node.seq.variable_indel.clear();
    node.seq.sequence = assign_sequence(node, alphabet);
    node.seq.non_char = range_union(&[node.seq.gaps.clone(), node.seq.unknown.clone()]);
    return Ok(vec![]);
  };

  if !is_leaf {
    node.seq.sequence = assign_sequence(node, alphabet);
  }
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

pub(crate) trait IndexedPartition {
  fn inputs(&self) -> &DenseInputs;

  fn indexed_kind(&self) -> IndexedKind<'_>;
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

impl IndexedKind<'_> {
  fn site_count(&self) -> usize {
    match self {
      IndexedKind::Dense { length, .. } => *length,
      IndexedKind::Discrete { .. } => 1,
    }
  }
}
