use crate::alphabet::alphabet::Alphabet;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::partition::marginal::dense::partition::assign_sequence;
use crate::partition::marginal::shared::data::{DenseInputs, count_transitions_dense};
use crate::partition::marginal::shared::normalize::{
  divide_out, forward_log_lh_add_normalization, forward_log_lh_remove_child, normalize_inplace,
  normalize_sites_from_log, normalize_sites_inplace,
};
use crate::partition::marginal::shared::update::{MarginalBackward, MarginalForward, MarginalPasses, TransitionCounts};
use crate::partition::storage::dense::{
  DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseLeafObs, DenseNodeState, DenseSeqDistribution,
  DenseSeqInfo,
};
use crate::seq::indel::{InDel, compute_node_ranges, resolve_indels_backward, resolve_indels_forward};
use eyre::Report;
use itertools::Itertools;
use ndarray::{Array1, Array2};
use parking_lot::Mutex;
use smallvec::smallvec;
use std::sync::Arc;
use treetime_graph::graph::Graph;
use treetime_graph::key_map::{EdgeMap, NodeMap};
use treetime_graph::node::GraphNodeKey;
use treetime_graph::schedule::GraphSchedule;
use treetime_graph::schedule_scan::{NodeEdgeMaps, PerEdge};
use treetime_graph::scope::Inheritance;
use treetime_primitives::LogLh;
use treetime_utils::interval::range_union::range_union;
use treetime_utils::{make_internal_error, make_internal_report};

impl<P: IndexedPartition> MarginalPasses for P {
  type Node = DenseNodeState;
  type BackwardInput = ();
  type Backward = DenseEdgeBackward;
  type Forward = DenseEdgeForward;
  type Estimate = DenseEdgeEstimate;

  fn inheritance(&self) -> &Inheritance {
    &self.inputs().inheritance
  }

  fn backward_input(_node_states: &NodeMap<DenseNodeState>) -> &Self::BackwardInput {
    &()
  }

  fn fresh_backward_input(_node_states: &NodeMap<DenseNodeState>) -> Self::BackwardInput {}

  fn marginal_backward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &EdgeMap<f64>,
    (): &(),
  ) -> Result<MarginalBackward<DenseNodeState, DenseEdgeBackward>, Report> {
    let kind = self.indexed_kind();
    indexed_backward(self.inputs(), gtr, kind, graph, branch_lengths)
  }

  fn marginal_forward(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &EdgeMap<f64>,
    node_states: NodeMap<DenseNodeState>,
    backward: &EdgeMap<DenseEdgeBackward>,
  ) -> Result<MarginalForward<DenseNodeState, DenseEdgeForward, DenseEdgeEstimate>, Report> {
    let kind = self.indexed_kind();
    indexed_forward(self.inputs(), gtr, kind, graph, branch_lengths, node_states, backward)
  }
}

impl<P: IndexedPartition> TransitionCounts for P {
  fn count_transitions(
    &self,
    gtr: &GTR,
    graph: &Graph,
    branch_lengths: &EdgeMap<f64>,
    node_states: &NodeMap<DenseNodeState>,
    backward: &EdgeMap<DenseEdgeBackward>,
    forward: &EdgeMap<DenseEdgeForward>,
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
  branch_lengths: &EdgeMap<f64>,
) -> Result<MarginalBackward<DenseNodeState, DenseEdgeBackward>, Report> {
  run_indexed_backward(inputs, gtr, kind, graph, branch_lengths, None)
}

pub(crate) fn indexed_backward_with_site_log_lh(
  partition: &impl IndexedPartition,
  gtr: &GTR,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  length: usize,
) -> Result<(MarginalBackward<DenseNodeState, DenseEdgeBackward>, Array1<f64>), Report> {
  let site_log_norms = Mutex::new(NodeMap::new());
  let backward = run_indexed_backward(
    partition.inputs(),
    gtr,
    partition.indexed_kind(),
    graph,
    branch_lengths,
    Some(&site_log_norms),
  )?;
  let site_log_lh = site_log_norms
    .into_inner()
    .into_values()
    .fold(Array1::zeros(length), |total, node| total + node);
  Ok((backward, site_log_lh))
}

fn run_indexed_backward(
  inputs: &DenseInputs,
  gtr: &GTR,
  kind: IndexedKind<'_>,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  site_log_norms: Option<&Mutex<NodeMap<Array1<f64>>>>,
) -> Result<MarginalBackward<DenseNodeState, DenseEdgeBackward>, Report> {
  let schedule = GraphSchedule::new(graph)?;
  let inheritance = &inputs.inheritance;
  inheritance.require_one_parent(&schedule)?;
  let pass = IndexedBackwardPass {
    schedule: &schedule,
    gtr,
    min_branch_length: inputs.min_branch_length,
    branch_lengths,
    site_log_norms,
    inheritance,
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
  let NodeEdgeMaps { nodes, edges } = inheritance.split_parent_edge(&schedule, outputs.unwrap_present().into_map()?)?;
  Ok(MarginalBackward {
    node_states: nodes,
    backward: edges,
  })
}

fn leaf_inputs<'a, L>(
  schedule: &'a GraphSchedule<'_>,
  leaves: &'a NodeMap<L>,
) -> impl Iterator<Item = (GraphNodeKey, Option<&'a L>)> {
  schedule.preorder().iter().map(|&key| (key, leaves.get(key)))
}

struct IndexedBackwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  gtr: &'a GTR,
  min_branch_length: f64,
  branch_lengths: &'a EdgeMap<f64>,
  site_log_norms: Option<&'a Mutex<NodeMap<Array1<f64>>>>,
  inheritance: &'a Inheritance,
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
    if !self.inheritance.holds(self.schedule, key) {
      return Ok((None, None));
    }
    let mut own_site_log_norms = None;
    let carried = self
      .inheritance
      .child_messages(self.schedule, key, children)?
      .into_iter()
      .map(|(_, _, child)| child.as_ref())
      .collect_vec();
    let (mut node, msg_to_parent) = if carried.is_empty() {
      let leaf =
        leaf.ok_or_else(|| make_internal_report!("Leaf {key} has no observation before the marginal backward pass"))?;
      leaf_backward(leaf)?
    } else {
      let child_nodes = carried.iter().map(|child| &child.0).collect_vec();
      let seq = internal_seq(&child_nodes);
      let child_edges = carried
        .iter()
        .map(|child| {
          child.1.first().ok_or_else(|| {
            make_internal_report!("A backward child of node {key} published no message for its parent edge")
          })
        })
        .collect::<Result<Vec<_>, Report>>()?;
      let Some((first_edge, other_edges)) = child_edges.split_first() else {
        return make_internal_error!("Internal node {key} has no children");
      };
      let mut log_dis = first_edge.msg_from_child.dis.mapv(f64::ln);
      for edge in other_edges {
        log_dis.zip_mut_with(&edge.msg_from_child.dis, |l, &p| *l += p.ln());
      }
      let (dis, site_log_norms) = normalize_sites_from_log(log_dis);
      let delta_ll = site_log_norms.iter().sum();
      own_site_log_norms = Some(site_log_norms);
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

    let parent_message = if let Some((_, _, edge_key)) = self.inheritance.parent(self.schedule, key) {
      let branch_length = self.branch_lengths[&edge_key].max(self.min_branch_length);
      smallvec![DenseEdgeBackward {
        msg_from_child: DenseSeqDistribution {
          dis: self.gtr.propagate_profile(&msg_to_parent.dis, branch_length, false),
          log_lh: msg_to_parent.log_lh,
        },
        msg_to_parent,
      }]
    } else {
      let mut dis = &msg_to_parent.dis * self.gtr.pi();
      let root_log_norms = normalize_sites_inplace(&mut dis);
      let delta_ll = root_log_norms.iter().sum();
      own_site_log_norms = Some(match own_site_log_norms {
        Some(site_log_norms) => site_log_norms + root_log_norms,
        None => root_log_norms,
      });
      node.profile = DenseSeqDistribution {
        dis,
        log_lh: msg_to_parent.log_lh + LogLh::new(delta_ll),
      };
      PerEdge::new()
    };

    if let (Some(site_log_norms), Some(own)) = (self.site_log_norms, own_site_log_norms) {
      site_log_norms.lock().insert(key, own);
    }
    let output = Arc::new((node, parent_message));
    Ok((Some(Arc::clone(&output)), Some(output)))
  }
}

type IndexedBackwardOutput = Option<Arc<(DenseNodeState, PerEdge<DenseEdgeBackward>)>>;

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
  branch_lengths: &EdgeMap<f64>,
  node_states: NodeMap<DenseNodeState>,
  backward: &EdgeMap<DenseEdgeBackward>,
) -> Result<MarginalForward<DenseNodeState, DenseEdgeForward, DenseEdgeEstimate>, Report> {
  let IndexedForwardMessages { node_states, forward } =
    indexed_forward_messages(inputs, gtr, graph, branch_lengths, node_states, backward)?;
  let (node_states, estimates) = indexed_forward_sequences(inputs, kind, graph, node_states)?;
  Ok(MarginalForward {
    node_states,
    forward,
    estimates,
  })
}

pub(crate) fn indexed_forward_messages(
  inputs: &DenseInputs,
  gtr: &GTR,
  graph: &Graph,
  branch_lengths: &EdgeMap<f64>,
  node_states: NodeMap<DenseNodeState>,
  backward: &EdgeMap<DenseEdgeBackward>,
) -> Result<IndexedForwardMessages, Report> {
  let schedule = GraphSchedule::new(graph)?;
  let inheritance = &inputs.inheritance;
  inheritance.require_one_parent(&schedule)?;
  let pass = IndexedForwardPass {
    schedule: &schedule,
    gtr,
    min_branch_length: inputs.min_branch_length,
    branch_lengths,
    inheritance,
  };
  let outputs = schedule.scan_down(
    inheritance.join_parent_edge(&schedule, node_states, backward.iter().collect()),
    |key, (node, backward), parents: PerEdge<Option<IndexedForwardOutput>>| {
      let Some(node) = node else {
        return Ok((None, None));
      };
      let parent = inheritance.parent_message(&schedule, key, &parents)?;
      let (output, message) = pass.node(key, node, backward, parent.map(|(_, _, parent)| &parent.0))?;
      Ok((Some(output), Some(message)))
    },
  )?;
  let NodeEdgeMaps { nodes, edges } = inheritance.split_parent_edge(&schedule, outputs.unwrap_present().into_map()?)?;
  Ok(IndexedForwardMessages {
    node_states: nodes,
    forward: edges,
  })
}

pub(crate) struct IndexedForwardMessages {
  pub(crate) node_states: NodeMap<DenseNodeState>,
  pub(crate) forward: EdgeMap<DenseEdgeForward>,
}

pub(crate) fn indexed_forward_sequences(
  inputs: &DenseInputs,
  kind: IndexedKind<'_>,
  graph: &Graph,
  node_states: NodeMap<DenseNodeState>,
) -> Result<(NodeMap<DenseNodeState>, EdgeMap<DenseEdgeEstimate>), Report> {
  let schedule = GraphSchedule::new(graph)?;
  let inheritance = &inputs.inheritance;
  inheritance.require_one_parent(&schedule)?;
  let mut node_states = node_states;
  let node_inputs = schedule
    .preorder()
    .iter()
    .map(|&key| (key, node_states.remove(key)))
    .collect_vec();
  let outputs = schedule.scan_down(
    node_inputs,
    |key, node: Option<DenseNodeState>, parents: PerEdge<Option<IndexedSequenceOutput>>| {
      let Some(mut node) = node else {
        return Ok((None, None));
      };
      let parent = inheritance
        .parent_message(&schedule, key, &parents)?
        .map(|(_, _, parent)| &parent.0);
      let indels = match kind {
        IndexedKind::Dense { alphabet, .. } => {
          let is_leaf = inheritance.children(&schedule, key).next().is_none();
          forward_post_dense(is_leaf, parent, &mut node, alphabet)?
        },
        IndexedKind::Discrete { .. } => vec![],
      };
      let edges = parent.map(|_| DenseEdgeEstimate { indels }).into_iter().collect();
      let output = Arc::new((node, edges));
      Ok((Some(Arc::clone(&output)), Some(output)))
    },
  )?;
  let NodeEdgeMaps { nodes, edges } = inheritance.split_parent_edge(&schedule, outputs.unwrap_present().into_map()?)?;
  Ok((nodes, edges))
}

type IndexedSequenceOutput = Arc<(DenseNodeState, PerEdge<DenseEdgeEstimate>)>;

struct IndexedForwardPass<'a> {
  schedule: &'a GraphSchedule<'a>,
  gtr: &'a GTR,
  min_branch_length: f64,
  branch_lengths: &'a EdgeMap<f64>,
  inheritance: &'a Inheritance,
}

impl IndexedForwardPass<'_> {
  fn node(
    &self,
    key: GraphNodeKey,
    mut node: DenseNodeState,
    backward: Option<&DenseEdgeBackward>,
    parent: Option<&DenseNodeState>,
  ) -> Result<(IndexedForwardOutput, IndexedForwardOutput), Report> {
    let msg_to_child = match (self.inheritance.parent(self.schedule, key), backward, parent) {
      (Some((_, _, edge_key)), Some(backward), Some(parent)) => {
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

    let edges = msg_to_child
      .map(|msg_to_child| DenseEdgeForward { msg_to_child })
      .into_iter()
      .collect();
    let output = Arc::new((node, edges));
    Ok((Arc::clone(&output), output))
  }
}

type IndexedForwardOutput = Arc<(DenseNodeState, PerEdge<DenseEdgeForward>)>;

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
  if is_leaf {
    for unknown in &node.seq.unknown {
      node.seq.sequence[unknown.0..unknown.1].fill(alphabet.unknown());
    }
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
    leaves: &'a NodeMap<DenseLeafObs>,
  },
  Discrete {
    leaves: &'a NodeMap<Array2<f64>>,
  },
}
