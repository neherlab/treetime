use crate::alphabet::alphabet::Alphabet;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::Weighted;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::io::fasta::FastaRecord;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::graph_dense::{DenseSeqDis, DenseSeqEdge, DenseSeqInfo, DenseSeqNode};
use crate::representation::log_lh::graph_log_lh;
use crate::representation::partition_marginal_dense::PartitionMarginalDense;
use crate::representation::seq::Seq;
use crate::utils::container::get_exactly_one;
use crate::utils::interval::range_intersection::range_intersection;
use crate::{make_report, seq};
use eyre::Report;
use itertools::Itertools;
use log::debug;
use ndarray::prelude::*;
use ndarray_stats::QuantileExt;
use parking_lot::RwLock;
use std::sync::Arc;

// TODO: move this into Alphabet
// Turn a profile into a sequence
fn prof2seq(profile: &DenseSeqDis, alphabet: &Alphabet) -> Seq {
  // iterate over profile.dis and take the argmax of each row, lookup in alphabet.char to extract the character
  let mut seq = seq! {};
  for row in profile.dis.rows() {
    let argmax = row.argmax().unwrap();
    seq.push(alphabet.char(argmax));
  }
  seq
}

fn assign_sequence(seq_info: &DenseSeqNode, alphabet: &Alphabet) -> Seq {
  let mut seq = prof2seq(&seq_info.profile, alphabet);
  for gap in &seq_info.seq.gaps {
    seq[gap.0..gap.1].fill(alphabet.gap());
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

fn attach_seqs_to_graph(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  aln: &[FastaRecord],
) -> Result<(), Report> {
  for leaf in graph.get_leaves() {
    let leaf_key = leaf.read_arc().key();
    let mut leaf = leaf.read_arc().payload().write_arc();

    let leaf_name = leaf.name.as_ref().ok_or_else(|| {
      make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
    })?.to_owned();

    let leaf_fasta = aln
      .iter()
      .find(|fasta| fasta.seq_name == leaf_name)
      .ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;

    *leaf = NodeAncestral {
      name: Some(leaf_name),
      desc: leaf_fasta.desc.clone(),
    };

    partitions.iter().try_for_each(|partition| -> Result<(), Report> {
      let mut partition = partition.write_arc();
      let alphabet = &partition.alphabet.clone(); // TODO: avoid clone

      partition
        .nodes
        .insert(leaf_key, DenseSeqNode::new(&leaf_fasta.seq, alphabet)?);

      Ok(())
    })?;
  }

  for edge in graph.get_edges() {
    let edge_key = edge.read_arc().key();
    partitions.iter().try_for_each(|partition| -> Result<(), Report> {
      let mut partition = partition.write_arc();
      partition.edges.insert(edge_key, DenseSeqEdge::default());
      Ok(())
    })?;
  }

  Ok(())
}

/// Backward pass calculates ingroup profiles
fn marginal_dense_backward(graph: &GraphAncestral, partitions: &[Arc<RwLock<PartitionMarginalDense>>]) {
  graph.par_iter_breadth_first_backward(|mut node| {
    run_marginal_dense_backward(partitions, &mut node).unwrap();
    GraphTraversalContinuation::Continue
  });
}

fn run_marginal_dense_backward(
  partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  node: &mut GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  for partition in partitions {
    let mut partition = partition.write_arc();
    let alphabet = &partition.alphabet.clone(); // TODO: avoid clone
    let length = partition.length;

    let msg_to_parent = if node.is_leaf {
      let seq_info = &partition.nodes[&node.key];
      DenseSeqDis {
        dis: alphabet.seq2prof(&seq_info.seq.sequence)?,
        log_lh: 0.0,
      }
    } else {
      let gaps = node
        .child_keys
        .iter()
        .map(|(child_key, _)| partition.nodes[child_key].seq.gaps.clone()) // TODO: avoid cloning
        .collect_vec();

      let gaps = range_intersection(&gaps);

      let seq = DenseSeqInfo {
        gaps,
        ..DenseSeqInfo::default()
      };

      let node_data = DenseSeqNode {
        seq,
        profile: DenseSeqDis::default(),
      };
      partition.nodes.insert(node.key, node_data);

      let msgs = node
        .child_keys
        .iter()
        .map(|(_, edge_key)| {
          partition.edges[edge_key].msg_from_child.dis.view().to_owned() //FIXME: avoid copy
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
        .map(|(_, edge_key)| partition.edges[edge_key].msg_from_child.log_lh)
        .sum::<f64>();
      DenseSeqDis {
        dis,
        log_lh: log_lh + delta_ll,
      }
    };

    if node.is_root {
      let mut dis = &msg_to_parent.dis * &partition.gtr.pi;
      let delta_ll = normalize_inplace(&mut dis);

      let seq_info = partition.nodes.get_mut(&node.key).unwrap();
      seq_info.profile.dis = dis;
      seq_info.profile.log_lh = msg_to_parent.log_lh + delta_ll;
    } else {
      // what was calculated above is what is sent to the parent. we also calculate the propagated message to the parent (we need it in the forward pass).
      let edge_key = get_exactly_one(&node.parent_edge_keys).expect("Only nodes with exactly one parent are supported");
      let branch_length = node.parent_edges[0].weight().unwrap_or(0.0);
      let branch_length = fix_branch_length(length, branch_length);
      let mut edge_data = DenseSeqEdge::default();

      let mut dis = Array2::ones((length, alphabet.n_canonical()));
      let log_lh = msg_to_parent.log_lh;
      let exp_qt = partition.gtr.expQt(branch_length);

      // Note that these dot products are really just
      //     a_{ik} = b.dot(c) = sum(b_{ij}c{jk}, j)
      // -- we can use whatever memory layout we want.
      dis *= &msg_to_parent.dis.dot(&exp_qt);
      // if let Some(transmission) = &edge_data.transmission {
      //   for r in transmission {
      //     dis
      //       .slice_mut(s![r.0..r.1, ..])
      //       .assign(&msg_to_parent.dis.slice(s![r.0..r.1, ..]).dot(&exp_qt));
      //   }
      // } else {
      //   // could make this a copy
      //   dis *= &msg_to_parent.dis.dot(&exp_qt);
      // }
      edge_data.msg_from_child = DenseSeqDis { dis, log_lh };
      edge_data.msg_to_parent = msg_to_parent;
      partition.edges.insert(*edge_key, edge_data);
    }
  }
  Ok(())
}

fn marginal_dense_forward(graph: &GraphAncestral, partitions: &[Arc<RwLock<PartitionMarginalDense>>]) {
  graph.par_iter_breadth_first_forward(|mut node| {
    run_marginal_dense_forward(graph, partitions, &mut node);
    GraphTraversalContinuation::Continue
  });
}

fn run_marginal_dense_forward(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  node: &mut GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
) {
  for partition in partitions {
    let mut partition = partition.write_arc();
    let length = partition.length;
    if !node.is_root {
      let mut seq_info = partition.nodes.remove(&node.key).unwrap();
      let mut msgs_to_combine: Vec<Array2<f64>> = vec![];
      let mut log_lh = 0.0;
      for (_, edge_key) in &node.parent_keys {
        let edge = &partition.edges[edge_key];
        let edge_payload = graph.get_edge(*edge_key).unwrap().read_arc().payload().read_arc();
        let exp_qt_matrix = partition.gtr.expQt(edge_payload.branch_length.unwrap_or(0.0));
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
      partition.nodes.insert(node.key, seq_info);
    }

    for child_edge_key in &node.child_edge_keys {
      let mut child_edge_data = partition.edges.remove(child_edge_key).unwrap();
      let node_data = &partition.nodes[&node.key];

      // this normalization isn't strictly necessary
      let mut dis = &node_data.profile.dis / &child_edge_data.msg_from_child.dis;
      let delta_ll = normalize_inplace(&mut dis);
      child_edge_data.msg_to_child = DenseSeqDis {
        dis,
        log_lh: node_data.profile.log_lh - child_edge_data.msg_from_child.log_lh + delta_ll,
      };
      partition.edges.insert(*child_edge_key, child_edge_data);
    }
  }
}

pub fn run_marginal_dense(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  aln: &[FastaRecord],
) -> Result<f64, Report> {
  attach_seqs_to_graph(graph, partitions, aln)?;
  marginal_dense_backward(graph, partitions);
  let log_lh = graph_log_lh(graph, partitions)?;
  debug!("Log likelihood: {log_lh}");
  marginal_dense_forward(graph, partitions);
  Ok(log_lh)
}

pub fn ancestral_reconstruction_marginal_dense(
  graph: &GraphAncestral,
  include_leaves: bool,
  partitions: &[Arc<RwLock<PartitionMarginalDense>>],
  mut visitor: impl FnMut(&NodeAncestral, &Seq),
) -> Result<(), Report> {
  graph.iter_depth_first_preorder_forward(|node| {
    if !include_leaves && node.is_leaf {
      return;
    }

    let seq = partitions
      .iter()
      .flat_map(|partition| {
        let partition = partition.read_arc();
        let seq_info = &partition.nodes[&node.key];
        if node.is_leaf {
          seq_info.seq.sequence.clone()
        } else {
          assign_sequence(seq_info, &partition.alphabet)
        }
      })
      .collect();

    visitor(&node.payload, &seq);
  });

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::commands::ancestral::fitch::get_common_length;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::json::{JsonPretty, json_write_str};
  use crate::io::nwk::nwk_read_str;
  use crate::representation::graph_ancestral::GraphAncestral;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  #[test]
  fn test_ancestral_reconstruction_marginal_dense() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
      indoc! {r#"
      >root
      TCAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >CD
      CCGGCCATGTATTG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#},
      &NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      TCGGCGCTGTATTGAC
      >AB
      ACATCGCTGTA-TGAC
      >CD
      TCGGCGGTGTATTG--
    "#},
      &NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let treat_gap_as_unknown = true;
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: true,
      ..JC69Params::default()
    })?;

    let partitions_marginal_dense = [Arc::new(RwLock::new(PartitionMarginalDense {
      index: 0,
      gtr,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    run_marginal_dense(&graph, &partitions_marginal_dense, &aln)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal_dense(&graph, false, &partitions_marginal_dense, |node, seq| {
      actual.insert(node.name.clone(), seq.to_string());
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  //   #[test]
  //   fn test_root_state() -> Result<(), Report> {
  //     rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
  //
  //     let aln = read_many_fasta_str(
  //       indoc! {r#"
  //       >root
  //       ACAGCCATGTATTG--
  //       >AB
  //       ACATCCCTGTA-TG--
  //       >A
  //       ACATCGCCNNA--GAC
  //       >B
  //       GCATCCCTGTA-NG--
  //       >CD
  //       CCGGCCATGTATTG--
  //       >C
  //       CCGGCGATGTRTTG--
  //       >D
  //       TCGGCCGTGTRTTG--
  //     "#},
  //       &NUC_ALPHABET,
  //     )?;
  //     let graph: DenseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
  //
  //     let treat_gap_as_unknown = true;
  //     let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
  //
  //     // use non-trivial GTR with non-uniform stationary distribution (tests correct use of transposed matrices)
  //     let mu = 1.0;
  //     let pi = array![0.2, 0.3, 0.15, 0.35];
  //     let gtr = GTR::new(GTRParams {
  //       alphabet: Alphabet::default(),
  //       W: None,
  //       pi,
  //       mu,
  //     })?;
  //
  //     let partitions = vec![PartitionLikelihoodWithAln::new(gtr, alphabet, aln)?];
  //
  //     let log_lh = run_marginal_dense(&graph, partitions, true)?;
  //     // from test_scripts/ancestral_dense.py
  //     pretty_assert_ulps_eq!(-59.20297892181229, log_lh, epsilon = 1e-6);
  //
  //     // test variable position distribution at the root for position 0 (from test_scripts/ancestral_dense.py)
  //     let pos_zero_root = array![0.28212327, 0.21643546, 0.13800802, 0.36343326];
  //     let root = &graph
  //       .get_exactly_one_root()?
  //       .read_arc()
  //       .payload()
  //       .read_arc()
  //       .dense_partitions[0];
  //     let pos: usize = 0;
  //     pretty_assert_ulps_eq!(root.profile.dis.slice(s![pos, 0..4]), &pos_zero_root, epsilon = 1e-6);
  //
  //     // pull out internal node AB for testing
  //     let node_ab = &graph
  //       .get_node(GraphNodeKey(1))
  //       .unwrap()
  //       .read_arc()
  //       .payload()
  //       .read_arc()
  //       .dense_partitions[0];
  //
  //     // test variable position distribution at internal node (from test_scripts/ancestral_dense.py)
  //     let pos: usize = 0;
  //     let pos_zero_ab = array![0.51275208, 0.09128506, 0.24647255, 0.14949031];
  //     pretty_assert_ulps_eq!(node_ab.profile.dis.slice(s![pos, 0..4]), &pos_zero_ab, epsilon = 1e-6);
  //
  //     // test variable position distribution at internal node (obtained from python treetime)
  //     let dis_ab = array![
  //       0.0013914677323952813,
  //       0.002087201598592933,
  //       0.042827146239885545,
  //       0.9536941844291262
  //     ];
  //     let pos: usize = 3;
  //     pretty_assert_ulps_eq!(node_ab.profile.dis.slice(s![pos, 0..4]), &dis_ab, epsilon = 1e-6);
  //
  //     // test whether the log likelihood is the same regardless of the root (here for node AB)
  //     let log_lh_ab = node_ab.profile.log_lh;
  //     pretty_assert_ulps_eq!(log_lh_ab, log_lh, epsilon = 1e-8);
  //
  //     Ok(())
  //   }
}
