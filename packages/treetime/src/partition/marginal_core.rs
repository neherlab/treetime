use crate::constants::SUPERTINY_NUMBER;
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::{
  MutationCounts, accumulate_mutation_counts, get_branch_mutation_matrix, is_profile_informative,
};
use crate::partition::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
use eyre::Report;
use itertools::izip;
use ndarray::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdge, GraphEdgeKey, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::array::softmax_with_log_norm::softmax_with_log_norm;
use treetime_utils::collections::container::get_exactly_one;

#[derive(Clone, Debug)]
pub struct MarginalData {
  pub gtr: GTR,
  pub nodes: BTreeMap<GraphNodeKey, DenseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, DenseEdgePartition>,
  pub min_branch_length: f64,
}

pub trait MarginalPartition<N, E>: Send + Sync
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn marginal_data(&self) -> &MarginalData;
  fn marginal_data_mut(&mut self) -> &mut MarginalData;

  fn leaf_profile(&self, node_key: GraphNodeKey) -> Result<DenseSeqDistribution, Report>;

  fn backward_internal_pre(&mut self, _node: &GraphNodeBackward<N, E, ()>) {}

  fn forward_post(&mut self, _graph: &Graph<N, E, ()>, _node: &GraphNodeForward<N, E, ()>) -> Result<(), Report> {
    Ok(())
  }
}

pub fn marginal_process_node_backward<N, E>(
  partition: &mut impl MarginalPartition<N, E>,
  node: &GraphNodeBackward<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let msg_to_parent = if node.is_leaf {
    partition.leaf_profile(node.key)?
  } else {
    partition.backward_internal_pre(node);

    let data = partition.marginal_data_mut();

    data.nodes.entry(node.key).or_insert_with(|| DenseNodePartition {
      seq: DenseSeqInfo::default(),
      profile: DenseSeqDistribution::default(),
    });

    let mut child_iter = node.child_keys.iter();
    let (_, first_edge_key) = child_iter.next().expect("Internal node must have children");
    let mut log_dis = data.edges[first_edge_key].msg_from_child.dis.mapv(f64::ln);
    for (_, edge_key) in child_iter {
      log_dis += &data.edges[edge_key].msg_from_child.dis.mapv(f64::ln);
    }

    let (dis, delta_ll) = normalize_from_log(&log_dis);
    let log_lh: f64 = node
      .child_keys
      .iter()
      .map(|(_, edge_key)| data.edges[edge_key].msg_from_child.log_lh)
      .sum();

    let node_data = data.nodes.get_mut(&node.key).unwrap();
    node_data.profile = DenseSeqDistribution {
      dis: dis.clone(),
      log_lh: log_lh + delta_ll,
    };

    DenseSeqDistribution {
      dis,
      log_lh: log_lh + delta_ll,
    }
  };

  let data = partition.marginal_data_mut();

  if node.is_root {
    let mut dis = &msg_to_parent.dis * &data.gtr.pi;
    let delta_ll = normalize_inplace(&mut dis);

    let node_data = data.nodes.get_mut(&node.key).unwrap();
    node_data.profile.dis = dis;
    node_data.profile.log_lh = msg_to_parent.log_lh + delta_ll;
  } else {
    let edge_key = get_exactly_one(&node.parent_edge_keys).expect("Only nodes with exactly one parent are supported");
    let branch_length = node.parent_edges[0].branch_length().unwrap_or(0.0);
    let branch_length = data.effective_branch_length(branch_length);

    let mut edge_data = data.edges.remove(edge_key).unwrap();

    let log_lh = msg_to_parent.log_lh;
    let dis = data.gtr.propagate_profile(&msg_to_parent.dis, branch_length, false);
    edge_data.msg_from_child = DenseSeqDistribution { dis, log_lh };
    edge_data.msg_to_parent = msg_to_parent;
    data.edges.insert(*edge_key, edge_data);
  }

  Ok(())
}

pub fn marginal_process_node_forward<N, E>(
  partition: &mut impl MarginalPartition<N, E>,
  graph: &Graph<N, E, ()>,
  node: &GraphNodeForward<N, E, ()>,
) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  if !node.is_root {
    let data = partition.marginal_data();
    let mut dis: Option<Array2<f64>> = None;
    let mut log_lh = 0.0;
    for (_, edge_key) in &node.parent_keys {
      let edge = &data.edges[edge_key];
      let edge_payload = graph.get_edge(*edge_key).unwrap().read_arc().payload().read_arc();
      let branch_length = edge_payload.branch_length().unwrap_or(0.0);
      let branch_length = data.effective_branch_length(branch_length);
      let msg_child = data.gtr.evolve(&edge.msg_to_child.dis, branch_length, false);
      log_lh += edge.msg_to_parent.log_lh + edge.msg_to_child.log_lh;
      dis = Some(match dis {
        None => &edge.msg_to_parent.dis * &msg_child,
        Some(mut d) => {
          d *= &edge.msg_to_parent.dis;
          d *= &msg_child;
          d
        },
      });
    }
    let mut dis = dis.expect("Non-root node must have at least one parent");
    let delta_ll = normalize_inplace(&mut dis);
    log_lh += delta_ll;

    let data = partition.marginal_data_mut();
    data.nodes.get_mut(&node.key).unwrap().profile = DenseSeqDistribution { dis, log_lh };
  }

  partition.forward_post(graph, node)?;

  let data = partition.marginal_data_mut();
  let node_profile_dis = data.nodes[&node.key].profile.dis.clone();
  let node_profile_log_lh = data.nodes[&node.key].profile.log_lh;

  for child_edge_key in &node.child_edge_keys {
    let child_edge_data = data.edges.get_mut(child_edge_key).unwrap();

    let safe_child = child_edge_data.msg_from_child.dis.mapv(|v| v.max(f64::MIN_POSITIVE));
    let mut dis = &node_profile_dis / &safe_child;
    let delta_ll = normalize_inplace(&mut dis);
    child_edge_data.msg_to_child = DenseSeqDistribution {
      dis,
      log_lh: node_profile_log_lh - child_edge_data.msg_from_child.log_lh + delta_ll,
    };
  }

  Ok(())
}

pub fn normalize_inplace(dis: &mut Array2<f64>) -> f64 {
  let norm = dis.sum_axis(Axis(1));
  let n_cols = dis.ncols() as f64;
  let mut log_lh = 0.0;
  for (ri, mut row) in dis.outer_iter_mut().enumerate() {
    let n = norm[ri];
    if n > 0.0 && n.is_finite() {
      row /= n;
      log_lh += n.ln();
    } else {
      row.fill(1.0 / n_cols);
      log_lh += f64::NEG_INFINITY;
    }
  }
  log_lh
}

impl MarginalData {
  pub fn effective_branch_length(&self, raw: f64) -> f64 {
    raw.max(self.min_branch_length)
  }

  /// Count posterior-weighted transitions from dense profile matrices.
  ///
  /// Shared by dense and discrete partitions (both store full profile matrices).
  pub fn count_transitions<N, E>(&self, graph: &Graph<N, E, ()>) -> Result<MutationCounts, Report>
  where
    N: GraphNode,
    E: GraphEdge + HasBranchLength,
  {
    let n_states = self.gtr.pi.len();
    let mut nij = Array2::zeros((n_states, n_states));
    let mut Ti = Array1::zeros(n_states);

    for edge in graph.get_edges() {
      let edge_arc = edge.read_arc();
      let branch_length = self.effective_branch_length(edge_arc.payload().read_arc().branch_length().unwrap_or(0.0));
      let edge_key = edge_arc.key();

      let edge_data = &self.edges[&edge_key];

      let exp_qt = self.gtr.expQt(branch_length) + SUPERTINY_NUMBER;
      let mut_stack = get_branch_mutation_matrix(&edge_data.msg_to_child.dis, &edge_data.msg_to_parent.dis, &exp_qt);
      accumulate_mutation_counts(&mut_stack, branch_length, &mut nij, &mut Ti);
    }

    let root = graph.get_exactly_one_root()?;
    let root_key = root.read_arc().key();
    let root_profile = &self.nodes[&root_key].profile.dis;
    let mut root_state = Array1::zeros(n_states);
    for row in root_profile.rows() {
      if is_profile_informative(row, n_states) {
        if let Some(root_idx) = argmax_first(&row) {
          root_state[root_idx] += 1.0;
        }
      }
    }

    nij.diag_mut().fill(0.0);

    Ok(MutationCounts { nij, Ti, root_state })
  }
}

pub fn normalize_from_log(log_dis: &Array2<f64>) -> (Array2<f64>, f64) {
  let mut dis = Array2::zeros(log_dis.raw_dim());
  let mut total_log_lh = 0.0;

  for (mut out_row, log_row) in izip!(dis.rows_mut(), log_dis.rows()) {
    let (normalized, log_norm) = softmax_with_log_norm(log_row);
    out_row.assign(&normalized);
    total_log_lh += log_norm;
  }

  (dis, total_log_lh)
}
