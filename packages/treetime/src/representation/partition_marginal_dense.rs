use super::graph_dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDis, DenseSeqInfo};
use crate::graph::edge::{GraphEdge, GraphEdgeKey, Weighted};
use crate::graph::graph::{Graph, GraphNodeBackward, GraphNodeForward};
use crate::graph::node::{Described, GraphNode, GraphNodeKey, Named};
use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::PartitionWithGtrInference;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::io::fasta::FastaRecord;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::log_lh::HasLogLh;
use crate::representation::partition_marginal::{PartitionMarginal, PartitionMarginalOps};
use crate::representation::seq::Seq;
use crate::seq::composition::Composition;
use crate::seq::mutation::Sub;
use crate::{alphabet::alphabet::Alphabet, make_report, seq};
use eyre::Report;
use itertools::Itertools;
use ndarray::prelude::*;
use ndarray_stats::QuantileExt;
use std::collections::BTreeMap;
use treetime_utils::container::get_exactly_one;
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
  pub fn get_sequence_length(&self) -> Option<usize> {
    Some(self.length)
  }
}

impl HasLogLh for PartitionMarginalDense {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }
}

impl PartitionMarginal for PartitionMarginalDense {}

impl<N, E> crate::commands::timetree::partition_ops::PartitionTimetreeOps<N, E> for PartitionMarginalDense
where
  N: GraphNode + Named,
  E: GraphEdge + Weighted,
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
  N: GraphNode + Named + Described,
  E: GraphEdge + Weighted,
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

  fn process_node_backward(
    &mut self,
    node: &GraphNodeBackward<N, E, ()>,
  ) -> Result<(), Report> {
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

      let msgs = node
        .child_keys
        .iter()
        .map(|(_, edge_key)| {
          self.edges[edge_key].msg_from_child.dis.view().to_owned() //FIXME: avoid copy
        })
        .collect_vec();

      let mut dis = msgs[0].clone().to_owned();
      for msg in &msgs[1..] {
        dis *= msg;
      }
      let delta_ll = normalize_inplace(&mut dis);
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
      let branch_length = node.parent_edges[0].weight().unwrap_or(0.0);
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

  fn process_node_forward(
    &mut self,
    graph: &Graph<N, E, ()>,
    node: &GraphNodeForward<N, E, ()>,
  ) -> Result<(), Report> {
    if !node.is_root {
      let mut seq_info = self.nodes.remove(&node.key).unwrap();
      let mut msgs_to_combine: Vec<Array2<f64>> = vec![];
      let mut log_lh = 0.0;
      for (_, edge_key) in &node.parent_keys {
        let edge = &self.edges[edge_key];
        let edge_payload = graph.get_edge(*edge_key).unwrap().read_arc().payload().read_arc();
        let branch_length = edge_payload.weight().unwrap_or(0.0);
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

  fn extract_ancestral_sequence(&mut self, node_key: GraphNodeKey) -> Seq {
    if let Some(seq_info) = self.nodes.get(&node_key) {
      assign_sequence(seq_info, &self.alphabet)
    } else {
      seq! {}
    }
  }

  fn reconstruct_node_sequence(
    &mut self,
    node: &GraphNodeForward<N, E, ()>,
    include_leaves: bool,
  ) -> Option<Seq> {
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

  fn get_sequence_length(&self) -> Option<usize> {
    Some(self.length)
  }
}

fn assign_sequence(seq_info: &DenseNodePartition, alphabet: &Alphabet) -> Seq {
  let mut seq = prof2seq(&seq_info.profile, alphabet);
  for gap in &seq_info.seq.gaps {
    seq[gap.0..gap.1].fill(alphabet.gap());
  }
  seq
}

fn prof2seq(profile: &DenseSeqDis, alphabet: &Alphabet) -> Seq {
  let mut seq = seq! {};
  for row in profile.dis.rows() {
    let argmax = row.argmax().unwrap();
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

impl PartitionWithGtrInference for PartitionMarginalDense {
  fn alphabet(&self) -> &Alphabet {
    &self.alphabet
  }

  fn get_seq_composition(&self, _node_key: GraphNodeKey) -> &Composition {
    unimplemented!("Dense node composition is not implemented yet")
  }

  fn get_edge_substitutions(&self, _edge_key: GraphEdgeKey, _graph: &GraphAncestral) -> Vec<Sub> {
    unimplemented!("Dense node substitutions lookup is not implemented yet")
  }
}
