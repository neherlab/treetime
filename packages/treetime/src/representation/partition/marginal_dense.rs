use crate::alphabet::alphabet::Alphabet;
use crate::commands::timetree::partition_ops::PartitionRerootOps;
use crate::gtr::gtr::GTR;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::make_report;
use crate::representation::partition::traits::HasLogLh;
use crate::representation::partition::traits::{PartitionMarginal, PartitionMarginalOps};
use crate::representation::payload::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDis, DenseSeqInfo};
use eyre::Report;
use itertools::Itertools;
use ndarray::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::edge::{EdgeOptimizeOps, GraphEdgeKey};
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, GraphNodeKey, Named, NodeAncestralOps};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{Seq, seq};
use treetime_utils::array::ndarray::argmax_first;
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range_intersection::range_intersection;

#[derive(Clone, Debug)]
pub struct PartitionMarginalDense {
  pub index: usize,
  pub gtr: GTR,
  pub alphabet: Alphabet,
  pub length: usize,
  pub nodes: BTreeMap<GraphNodeKey, DenseNodePartition>,
  pub edges: BTreeMap<GraphEdgeKey, DenseEdgePartition>,
}

impl PartitionMarginalDense {
  #[allow(clippy::same_name_method)]
  pub fn get_sequence_length(&self) -> usize {
    self.length
  }
}

impl HasLogLh for PartitionMarginalDense {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }
}

impl PartitionMarginal for PartitionMarginalDense {}

impl PartitionRerootOps for PartitionMarginalDense {}

impl<N, E> crate::commands::timetree::partition_ops::PartitionTimetreeOps<N, E> for PartitionMarginalDense
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn create_edge_contribution(
    &self,
    edge_key: GraphEdgeKey,
  ) -> Result<crate::commands::optimize::optimize_unified::OptimizationContribution, Report> {
    Ok(crate::commands::optimize::optimize_unified::OptimizationContribution::from_dense(edge_key, self))
  }
}

impl<N, E> PartitionMarginalOps<N, E> for PartitionMarginalDense
where
  N: NodeAncestralOps,
  E: EdgeOptimizeOps,
{
  fn attach_sequences(&mut self, graph: &Graph<N, E, ()>, aln: &[FastaRecord]) -> Result<(), Report> {
    for leaf in graph.get_leaves() {
      let leaf_key = leaf.read_arc().key();
      let mut leaf = leaf.read_arc().payload().write_arc();

      let leaf_name = {
        let name = leaf.name().ok_or_else(|| {
          make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
        })?;
        name.as_ref().to_owned()
      };

      let leaf_fasta = aln
        .iter()
        .find(|fasta| fasta.seq_name == leaf_name)
        .ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;

      leaf.set_desc(leaf_fasta.desc.clone());

      let alphabet = &self.alphabet.clone(); // TODO: avoid clone
      self
        .nodes
        .insert(leaf_key, DenseNodePartition::new(&leaf_fasta.seq, alphabet)?);
    }

    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      self.edges.insert(edge_key, DenseEdgePartition::default());
    }

    Ok(())
  }

  fn process_node_backward(&mut self, node: &GraphNodeBackward<N, E, ()>) -> Result<(), Report> {
    let alphabet = &self.alphabet.clone(); // TODO: avoid clone
    let length = self.length;

    let msg_to_parent = if node.is_leaf {
      let seq_info = &self.nodes[&node.key];
      DenseSeqDis {
        dis: alphabet.seq2prof(&seq_info.seq.sequence)?,
        log_lh: 0.0,
      }
    } else {
      let gaps = node
        .child_keys
        .iter()
        .map(|(child_key, _)| self.nodes[child_key].seq.gaps.clone()) // TODO: avoid cloning
        .collect_vec();

      let gaps = range_intersection(&gaps);

      let seq = DenseSeqInfo {
        gaps,
        ..DenseSeqInfo::default()
      };

      let node_data = DenseNodePartition {
        seq,
        profile: DenseSeqDis::default(),
      };
      self.nodes.insert(node.key, node_data);

      // Combine child messages in log-space to match v0's numerical behavior.
      // v0 uses: tmp_log_subtree_LH += ch.marginal_log_Lx (treeanc.py:874)
      // where marginal_log_Lx = np.log(profile.dot(Qt))
      let msgs = node
        .child_keys
        .iter()
        .map(|(_, edge_key)| {
          self.edges[edge_key].msg_from_child.dis.view().to_owned() //FIXME: avoid copy
        })
        .collect_vec();

      // Convert to log and sum (equivalent to multiply in prob space)
      let mut log_dis = msgs[0].mapv(f64::ln);
      for msg in &msgs[1..] {
        log_dis += &msg.mapv(f64::ln);
      }

      // Normalize using v0's normalize_profile(log=True) algorithm:
      // 1. Subtract row max (logsumexp trick)
      // 2. Exponentiate
      // 3. Normalize by row sum
      let (dis, delta_ll) = normalize_from_log(&log_dis);
      let log_lh = node
        .child_keys
        .iter()
        .map(|(_, edge_key)| self.edges[edge_key].msg_from_child.log_lh)
        .sum::<f64>();
      DenseSeqDis {
        dis,
        log_lh: log_lh + delta_ll,
      }
    };

    if node.is_root {
      let mut dis = &msg_to_parent.dis * &self.gtr.pi;
      let delta_ll = normalize_inplace(&mut dis);

      let seq_info = self.nodes.get_mut(&node.key).unwrap();
      seq_info.profile.dis = dis;
      seq_info.profile.log_lh = msg_to_parent.log_lh + delta_ll;
    } else {
      let edge_key = get_exactly_one(&node.parent_edge_keys).expect("Only nodes with exactly one parent are supported");
      let branch_length = node.parent_edges[0].branch_length().unwrap_or(0.0);
      let branch_length = fix_branch_length(length, branch_length);
      let mut edge_data = DenseEdgePartition::default();

      let mut dis = Array2::ones((length, alphabet.n_canonical()));
      let log_lh = msg_to_parent.log_lh;
      let exp_qt = self.gtr.expQt(branch_length);

      dis *= &msg_to_parent.dis.dot(&exp_qt);
      edge_data.msg_from_child = DenseSeqDis { dis, log_lh };
      edge_data.msg_to_parent = msg_to_parent;
      self.edges.insert(*edge_key, edge_data);
    }
    Ok(())
  }

  fn process_node_forward(&mut self, graph: &Graph<N, E, ()>, node: &GraphNodeForward<N, E, ()>) -> Result<(), Report> {
    if !node.is_root {
      let mut seq_info = self.nodes.remove(&node.key).unwrap();
      let mut msgs_to_combine: Vec<Array2<f64>> = vec![];
      let mut log_lh = 0.0;
      for (_, edge_key) in &node.parent_keys {
        let edge = &self.edges[edge_key];
        let edge_payload = graph.get_edge(*edge_key).unwrap().read_arc().payload().read_arc();
        let branch_length = edge_payload.branch_length().unwrap_or(0.0);
        let exp_qt_matrix = self.gtr.expQt(branch_length);
        let exp_qt = exp_qt_matrix.t();
        msgs_to_combine.push(edge.msg_to_parent.dis.view().to_owned()); // FIXME: avoid copy
        log_lh += edge.msg_to_parent.log_lh;
        msgs_to_combine.push(edge.msg_to_child.dis.dot(&exp_qt));
        log_lh += edge.msg_to_child.log_lh;
      }
      let mut dis = msgs_to_combine[0].clone();
      for msg in &msgs_to_combine[1..] {
        dis *= msg;
      }
      let delta_ll = normalize_inplace(&mut dis);
      log_lh += delta_ll;
      seq_info.profile = DenseSeqDis { dis, log_lh };
      self.nodes.insert(node.key, seq_info);
    }

    for child_edge_key in &node.child_edge_keys {
      let mut child_edge_data = self.edges.remove(child_edge_key).unwrap();
      let node_data = &self.nodes[&node.key];

      // this normalization isn't strictly necessary
      let mut dis = &node_data.profile.dis / &child_edge_data.msg_from_child.dis;
      let delta_ll = normalize_inplace(&mut dis);
      child_edge_data.msg_to_child = DenseSeqDis {
        dis,
        log_lh: node_data.profile.log_lh - child_edge_data.msg_from_child.log_lh + delta_ll,
      };
      self.edges.insert(*child_edge_key, child_edge_data);
    }
    Ok(())
  }

  fn extract_ancestral_sequence(&self, node_key: GraphNodeKey) -> Seq {
    if let Some(seq_info) = self.nodes.get(&node_key) {
      assign_sequence(seq_info, &self.alphabet)
    } else {
      seq! {}
    }
  }

  fn reconstruct_node_sequence(&mut self, node: &GraphNodeForward<N, E, ()>, include_leaves: bool) -> Option<Seq> {
    if !include_leaves && node.is_leaf {
      return None;
    }

    let seq_info = self.nodes.get(&node.key)?;
    let seq = if node.is_leaf {
      seq_info.seq.sequence.clone()
    } else {
      assign_sequence(seq_info, &self.alphabet)
    };

    Some(seq)
  }

  fn get_sequence_length(&self) -> usize {
    self.length
  }
}

fn assign_sequence(seq_info: &DenseNodePartition, alphabet: &Alphabet) -> Seq {
  let mut seq = prof2seq(&seq_info.profile, alphabet);
  for gap in &seq_info.seq.gaps {
    seq[gap.0..gap.1].fill(alphabet.gap());
  }
  seq
}

/// Extract sequence from profile by taking argmax at each position.
///
/// Uses `argmax_first` for NumPy-compatible tie-breaking: when multiple states
/// have equal probability, returns the first (lowest index) state.
fn prof2seq(profile: &DenseSeqDis, alphabet: &Alphabet) -> Seq {
  let mut seq = seq! {};
  for row in profile.dis.rows() {
    let argmax = argmax_first(&row).unwrap_or(0);
    seq.push(alphabet.char(argmax));
  }
  seq
}

fn normalize_inplace(dis: &mut Array2<f64>) -> f64 {
  let norm = dis.sum_axis(Axis(1));
  for (ri, mut row) in dis.outer_iter_mut().enumerate() {
    row /= norm[ri];
  }
  norm.mapv(f64::ln).sum()
}

/// Normalize a log-probability matrix to probability matrix.
///
/// Matches v0's normalize_profile(in_profile, log=True) algorithm (seq_utils.py:296-307):
/// 1. Subtract row max (logsumexp trick for numerical stability)
/// 2. Exponentiate
/// 3. Normalize by row sum
fn normalize_from_log(log_dis: &Array2<f64>) -> (Array2<f64>, f64) {
  let (nrows, ncols) = log_dis.dim();
  let mut dis = Array2::zeros((nrows, ncols));
  let mut total_log_lh = 0.0;

  for (ri, log_row) in log_dis.outer_iter().enumerate() {
    // Step 1: Find row max (logsumexp trick)
    let max_val = log_row.iter().fold(f64::NEG_INFINITY, |a, &b| a.max(b));

    // Step 2: Exponentiate (relative to max)
    let mut row_sum = 0.0;
    for (ci, &log_val) in log_row.iter().enumerate() {
      let val = (log_val - max_val).exp();
      dis[[ri, ci]] = val;
      row_sum += val;
    }

    // Step 3: Normalize by row sum
    for ci in 0..ncols {
      dis[[ri, ci]] /= row_sum;
    }

    // Accumulate log-likelihood: log(sum(exp(log_vals))) = max_val + log(row_sum)
    total_log_lh += max_val + row_sum.ln();
  }

  (dis, total_log_lh)
}
