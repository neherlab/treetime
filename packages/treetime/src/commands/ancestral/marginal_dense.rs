use crate::alphabet::alphabet::Alphabet;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::Weighted;
use crate::hacks::fix_branch_length::fix_branch_length;
use crate::representation::graph_dense::{
  DenseGraph, DenseNode, DenseSeqDis, DenseSeqEdge, DenseSeqInfo, DenseSeqNode,
};
use crate::representation::partitions_likelihood::{PartitionLikelihood, PartitionLikelihoodWithAln};
use crate::utils::container::get_exactly_one_mut;
use crate::utils::interval::range_intersection::range_intersection;
use crate::{make_internal_report, make_report};
use eyre::Report;
use itertools::Itertools;
use log::debug;
use ndarray::prelude::*;
use ndarray::AssignElem;
use ndarray_stats::QuantileExt;

// Turn a profile into a sequence
fn prof2seq(profile: &DenseSeqDis, alphabet: &Alphabet) -> Vec<char> {
  // iterate over profile.dis and take the argmax of each row, lookup in alphabet.char to extract the character
  let mut seq = Vec::new();
  for row in profile.dis.rows() {
    let argmax = row.argmax().unwrap();
    seq.push(alphabet.char(argmax));
  }
  seq
}

fn assign_sequence(seq_info: &DenseSeqNode, alphabet: &Alphabet) -> Vec<char> {
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
  norm.mapv(|x| x.ln()).sum()
}

fn attach_seqs_to_graph(graph: &DenseGraph, partitions: &[PartitionLikelihoodWithAln]) -> Result<(), Report> {
  graph.get_leaves().iter().try_for_each(|leaf| -> Result<(), Report>  {
    let mut leaf = leaf.write_arc().payload().write_arc();

    let leaf_name = leaf.name.as_ref().ok_or_else(|| {
      make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
    })?.to_owned();

    // FIXME: all descs are the same for fasta partitions, so the mutable assignment here is needlessly complicated
    let mut desc = None;

    leaf.dense_partitions = partitions
      .iter()
      .map(|PartitionLikelihoodWithAln { gtr, alphabet, aln, length }| {
        // TODO(perf): this might be slow if there are many sequences
        let leaf_fasta = aln
          .iter()
          .find(|fasta| fasta.seq_name == leaf_name)
          .ok_or_else(|| make_internal_report!("Leaf sequence not found: '{leaf_name}'"))?;

        desc.assign_elem(leaf_fasta.desc.clone());

        // TODO(perf): unnecessary copy of sequence data. Neither String, nor &[char] works well for us, it seems.
        // We probably want a custom class for sequences. Sequences should be instantiated in the fasta parser and
        // never need a copy like here.
        let sequence = leaf_fasta.seq.chars().collect::<Vec<_>>();

        DenseSeqNode::new(&sequence, alphabet)
      })
      .collect::<Result<_, Report>>()?;

    leaf.desc = desc;

    Ok(())
  })?;

  graph.get_edges().iter().for_each(|edge| {
    let mut edge = edge.write_arc().payload().write_arc();
    edge.dense_partitions = vec![];
  });

  graph.get_internal_nodes().iter().for_each(|node| {
    let mut node = node.write_arc().payload().write_arc();
    node.dense_partitions = vec![];
  });

  Ok(())
}

fn ingroup_profiles_dense(graph: &DenseGraph, partitions: &[PartitionLikelihood]) {
  let n_partitions = partitions.len();
  graph.par_iter_breadth_first_backward(|mut node| {
    for si in 0..n_partitions {
      let PartitionLikelihood { gtr, alphabet, length } = &partitions[si];
      let msg_to_parent = if node.is_leaf {
        let seq_info = &node.payload.dense_partitions[si];
        DenseSeqDis {
          dis: alphabet.seq2prof(&seq_info.seq.sequence).unwrap(),
          log_lh: 0.0,
        }
      } else {
        let gaps = node
        .children
        .iter()
        .map(|(c, e)| c.read_arc().dense_partitions[si].seq.gaps.clone()) // TODO: avoid cloning
        .collect_vec();

        let gaps = range_intersection(&gaps);

        let seq = DenseSeqInfo {
          gaps,
          ..DenseSeqInfo::default()
        };
        node.payload.dense_partitions.push(DenseSeqNode {
          seq,
          profile: DenseSeqDis::default(),
        });

        let msgs = node
          .children
          .iter()
          .map(|(c, e)| {
            let edge = e.read_arc();
            edge.dense_partitions[si].msg_from_child.dis.view().to_owned() //FIXME: avoid copy
          })
          .collect_vec();

        let mut dis = msgs[0].clone().to_owned();
        for msg in msgs[1..].iter() {
          dis *= msg;
        }
        let delta_ll = normalize_inplace(&mut dis);
        let log_lh = node
          .children
          .iter()
          .map(|(c, e)| {
            let edge = e.read_arc();
            edge.dense_partitions[si].msg_from_child.log_lh
          })
          .sum::<f64>();
        DenseSeqDis {
          dis,
          log_lh: log_lh + delta_ll,
        }
      };

      if node.is_root {
        let seq_info = &mut node.payload.dense_partitions[si];
        let mut dis = &msg_to_parent.dis * &gtr.pi;
        let delta_ll = normalize_inplace(&mut dis);

        seq_info.profile.dis = dis;
        seq_info.profile.log_lh = msg_to_parent.log_lh + delta_ll;
      } else {
        // what was calculated above is what is sent to the parent. we also calculate the propagated message to the parent (we need it in the forward pass).
        let edge_to_parent =
          get_exactly_one_mut(&mut node.parent_edges).expect("Only nodes with exactly one parent are supported"); // HACK
        let branch_length = edge_to_parent.weight().unwrap_or(0.0);
        let branch_length = fix_branch_length(*length, branch_length);
        let mut edge_data = DenseSeqEdge::default();

        let mut dis = Array2::ones((*length, alphabet.n_canonical()));
        let log_lh = msg_to_parent.log_lh;
        let exp_qt = gtr.expQt(branch_length);

        // Note that these dot products are really just
        //     a_{ik} = b.dot(c) = sum(b_{ij}c{jk}, j)
        // -- we can use whatever memory layout we want.
        dis *= &msg_to_parent.dis.dot(&exp_qt);
        // if let Some(transmission) = &edge_to_parent.dense_partitions[si].transmission {
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
        edge_to_parent.dense_partitions.push(edge_data);
      };
    }
    GraphTraversalContinuation::Continue
  });
}

fn outgroup_profiles_dense(graph: &DenseGraph, partitions: &[PartitionLikelihood]) {
  let n_partitions = partitions.len();
  graph.par_iter_breadth_first_forward(|mut node| {
    for si in 0..n_partitions {
      let PartitionLikelihood { gtr, alphabet, length } = &partitions[si];
      if !node.is_root {
        let seq_info = &mut node.payload.dense_partitions[si];
        let mut msgs_to_combine: Vec<Array2<f64>> = vec![];
        let mut log_lh = 0.0;
        for (pi, (p, edge)) in node.parents.iter().enumerate() {
          let edge = edge.read_arc();
          let expQt_matrix = gtr.expQt(edge.branch_length.unwrap_or(0.0));
          let expQt = expQt_matrix.t();
          msgs_to_combine.push(edge.dense_partitions[si].msg_to_parent.dis.view().to_owned()); //FIXME: avoid copy
          log_lh += edge.dense_partitions[si].msg_to_parent.log_lh;
          msgs_to_combine.push(edge.dense_partitions[si].msg_to_child.dis.dot(&expQt));
          log_lh += edge.dense_partitions[si].msg_to_child.log_lh;
        }
        let mut dis = msgs_to_combine[0].clone();
        for msg in msgs_to_combine[1..].iter() {
          dis *= msg;
        }
        let delta_ll = normalize_inplace(&mut dis);
        log_lh += delta_ll;
        seq_info.profile = DenseSeqDis { dis, log_lh };
      };

      for child_edge in &mut node.child_edges {
        // this normalization isn't strictly necessary
        let mut dis =
          &node.payload.dense_partitions[si].profile.dis / &child_edge.dense_partitions[si].msg_from_child.dis;
        let delta_ll = normalize_inplace(&mut dis);
        child_edge.dense_partitions[si].msg_to_child = DenseSeqDis {
          dis,
          log_lh: node.payload.dense_partitions[si].profile.log_lh
            - child_edge.dense_partitions[si].msg_from_child.log_lh
            + delta_ll,
        };
      }
    }
    GraphTraversalContinuation::Continue
  });
}

pub fn run_marginal_dense(
  graph: &DenseGraph,
  partitions: Vec<PartitionLikelihoodWithAln>,
  reconstruct: bool,
) -> Result<f64, Report> {
  attach_seqs_to_graph(graph, &partitions)?;

  let partitions = partitions.into_iter().map(PartitionLikelihood::from).collect_vec();

  ingroup_profiles_dense(graph, &partitions);
  let log_lh = graph
    .get_exactly_one_root()
    .unwrap()
    .read_arc()
    .payload()
    .read_arc()
    .dense_partitions
    .iter()
    .map(|p| p.profile.log_lh)
    .sum();
  debug!("Log likelihood: {}", log_lh);
  outgroup_profiles_dense(graph, &partitions);

  if reconstruct {
    for node in graph.get_nodes() {
      if node.read_arc().is_leaf() {
        continue;
      }
      for (si, seq_info) in node
        .write_arc()
        .payload()
        .write_arc()
        .dense_partitions
        .iter_mut()
        .enumerate()
      {
        let alphabet = &partitions[si].alphabet;
        seq_info.seq.sequence = assign_sequence(seq_info, alphabet);
      }
    }
  }
  Ok(log_lh)
}

pub fn ancestral_reconstruction_marginal_dense(
  graph: &DenseGraph,
  include_leaves: bool,
  mut visitor: impl FnMut(&DenseNode, Vec<char>),
) -> Result<(), Report> {
  graph.iter_depth_first_preorder_forward(|node| {
    if !include_leaves && node.is_leaf {
      return;
    }

    let seq = node
      .payload
      .dense_partitions
      .iter()
      .flat_map(|p| p.seq.sequence.iter().copied())
      .collect();

    visitor(&node.payload, seq);
  });

  Ok(())
}

#[cfg(test)]
mod tests {
  use std::collections::BTreeMap;

  use super::*;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::graph::node::GraphNodeKey;
  use crate::gtr::get_gtr::{jc69, JC69Params};
  use crate::gtr::gtr::{GTRParams, GTR};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::json::{json_write_str, JsonPretty};
  use crate::io::nwk::nwk_read_str;
  use crate::pretty_assert_ulps_eq;
  use crate::utils::string::vec_to_string;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use pretty_assertions::assert_eq;

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

    let graph: DenseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let treat_gap_as_unknown = true;
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;
    let gtr = jc69(JC69Params {
      alphabet: AlphabetName::Nuc,
      treat_gap_as_unknown: true,
      ..JC69Params::default()
    })?;
    let partitions = vec![PartitionLikelihoodWithAln::new(gtr, alphabet, aln)?];
    run_marginal_dense(&graph, partitions, true)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_marginal_dense(&graph, false, |node, seq| {
      actual.insert(node.name.clone(), vec_to_string(seq));
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[test]
  fn test_root_state() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
      indoc! {r#"
      >root
      ACAGCCATGTATTG--
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
    let graph: DenseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let treat_gap_as_unknown = true;
    let alphabet = Alphabet::new(AlphabetName::Nuc, treat_gap_as_unknown)?;

    let pi = array![0.2, 0.3, 0.15, 0.35];

    // symmetric rate matrix
    let mu = 1.0;

    let gtr = GTR::new(GTRParams {
      alphabet: Alphabet::default(),
      W: None,
      pi: pi.clone(),
      mu,
    })
    .unwrap();

    let partitions = vec![PartitionLikelihoodWithAln::new(gtr, alphabet, aln)?];

    let log_lh = run_marginal_dense(&graph, partitions, true)?;

    // test variable position distribution at the root
    let pos_zero_root = array![0.28212327, 0.21643546, 0.13800802, 0.36343326];
    let root = &graph
      .get_exactly_one_root()
      .unwrap()
      .read_arc()
      .payload()
      .read_arc()
      .dense_partitions[0];
    let pos: usize = 0;
    dbg!(&root.profile.dis.slice(s![pos, 0..4]));

    dbg!(log_lh);
    // pretty_assert_ulps_eq!(-56.946298878390444, log_lh, epsilon = 1e-6);
    pretty_assert_ulps_eq!(root.profile.dis.slice(s![pos, 0..4]), &pos_zero_root, epsilon = 1e-6);

    // test variable position distribution at internal node
    let pos_zero_ab = array![0.51275208, 0.09128506, 0.24647255, 0.14949031];
    let node_ab = &graph
      .get_node(GraphNodeKey(1))
      .unwrap()
      .read_arc()
      .payload()
      .read_arc()
      .dense_partitions[0];
    let pos: usize = 0;
    pretty_assert_ulps_eq!(node_ab.profile.dis.slice(s![pos, 0..4]), &pos_zero_ab, epsilon = 1e-6);

    // test variable position distribution at internal node
    let dis_ab = array![
      0.0013914677323952813,
      0.002087201598592933,
      0.042827146239885545,
      0.9536941844291262
    ];
    let pos: usize = 3;
    pretty_assert_ulps_eq!(node_ab.profile.dis.slice(s![pos, 0..4]), &dis_ab, epsilon = 1e-6);

    let log_lh_ab = node_ab.profile.log_lh;
    pretty_assert_ulps_eq!(log_lh_ab, log_lh, epsilon = 1e-8);

    Ok(())
  }
}
