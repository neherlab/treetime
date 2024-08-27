use crate::alphabet::alphabet::Alphabet;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::edge::Weighted;
use crate::graph::node::Named;
use crate::port::seq_dense::{DenseGraph, DenseNode, DenseSeqDis, DenseSeqEdge, DenseSeqInfo, DenseSeqNode};
use crate::port::seq_partitions::{PartitionLikelihood, PartitionLikelihoodWithAln};
use crate::seq::range_intersection::range_intersection;
use crate::utils::ndarray::{log, product_axis};
use crate::{make_internal_report, make_report, o};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use ndarray::prelude::*;
use ndarray::stack;
use ndarray_stats::QuantileExt;
use std::collections::BTreeMap;

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

fn attach_seqs_to_graph(graph: &DenseGraph, partitions: &[PartitionLikelihoodWithAln]) -> Result<(), Report> {
  graph.get_leaves().iter().try_for_each(|leaf| -> Result<(), Report>  {
    let mut leaf = leaf.write_arc().payload().write_arc();

    let leaf_name = leaf.name.as_ref().ok_or_else(|| {
      make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
    })?.to_owned();

    leaf.dense_partitions = partitions
      .iter()
      .map(|PartitionLikelihoodWithAln { gtr, alphabet, aln, length }| {
        // TODO(perf): this might be slow if there are many sequences
        let leaf_fasta = aln
          .iter()
          .find(|fasta| fasta.seq_name == leaf_name)
          .ok_or_else(|| make_internal_report!("Leaf sequence not found: '{leaf_name}'"))?;

        // TODO(perf): unnecessary copy of sequence data. Neither String, nor &[char] works well for us, it seems.
        // We probably want a custom class for sequences. Sequences should be instantiated in the fasta parser and
        // never need a copy like here.
        let sequence = leaf_fasta.seq.chars().collect::<Vec<_>>();

        DenseSeqNode::new(&sequence, alphabet)
      })
      .collect::<Result<_, Report>>()?;

    Ok(())
  })?;

  graph.get_edges().iter().for_each(|edge| {
    let mut edge = edge.write_arc().payload().write_arc();
    edge.dense_partitions = partitions.iter().map(|_| DenseSeqEdge::default()).collect_vec();
  });

  Ok(())
}

fn combine_dense_messages(msgs: &BTreeMap<String, DenseSeqDis>) -> Result<DenseSeqDis, Report> {
  let stacked = stack(Axis(0), &msgs.values().map(|d| d.dis.view()).collect_vec())?;
  let prod_dis = product_axis(&stacked, Axis(0));
  let norm = prod_dis.sum_axis(Axis(1));
  let dis = (&prod_dis.t() / &norm).t().to_owned();
  let log_lh = msgs.values().map(|m| m.log_lh).sum::<f64>() + log(&norm).sum();
  Ok(DenseSeqDis { dis, log_lh })
}

fn ingroup_profiles_dense(graph: &DenseGraph, partitions: &[PartitionLikelihood]) {
  let n_partitions = partitions.len();

  graph.par_iter_breadth_first_backward(|mut node| {
    if node.is_leaf {
      return GraphTraversalContinuation::Continue; // the msg to parents for leaves was put in init phase
    }

    node.payload.dense_partitions = (0..n_partitions)
      .map(|si| {
        let PartitionLikelihood { gtr, alphabet, length } = &partitions[si];

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

        let mut msgs_from_children = btreemap! {};

        for (child, edge) in &node.children {
          let (child, edge) = (child.read_arc(), edge.read_arc());

          let name = child
            .name()
            .expect("Encountered child node without a name")
            .as_ref()
            .to_owned();

          let exp_qt = gtr.expQt(edge.weight().unwrap_or_default());
          let child = &child.dense_partitions[si];
          let edge = &edge.dense_partitions[si];

          let mut dis = Array2::ones((*length, alphabet.n_canonical()));
          let log_lh = child.msg_to_parents.log_lh;

          // Note that these dot products are really just
          //     a_{ik} = b.dot(c) = sum(b_{ij}c{jk}, j)
          // -- we can use whatever memory layout we want.
          if let Some(transmission) = &edge.transmission {
            for r in transmission {
              dis
                .slice_mut(s![r.0..r.1, ..])
                .assign(&child.msg_to_parents.dis.slice(s![r.0..r.1, ..]).dot(&exp_qt));
            }
          } else {
            // could make this a copy
            dis *= &child.msg_to_parents.dis.dot(&exp_qt);
          }

          msgs_from_children.insert(name, DenseSeqDis { dis, log_lh });
        }

        let msg_to_parents = combine_dense_messages(&msgs_from_children).unwrap();

        DenseSeqNode {
          seq,
          msg_to_parents,
          msgs_from_children,
          ..DenseSeqNode::default()
        }
      })
      .collect_vec();

    GraphTraversalContinuation::Continue
  });
}

fn outgroup_profiles_dense(graph: &DenseGraph, dense_partitions: &[PartitionLikelihood]) {
  graph.par_iter_breadth_first_forward(|mut node| {
    let name = node
      .payload
      .name()
      .expect("Encountered node without a name")
      .as_ref()
      .to_owned();

    if node.is_root {
      return GraphTraversalContinuation::Continue;
    }

    for (si, seq_info) in node.payload.dense_partitions.iter_mut().enumerate() {
      let PartitionLikelihood { gtr, alphabet, length } = &dense_partitions[si];

      let mut msgs_from_parents = btreemap! {};
      for (parent, edge) in &node.parents {
        let parent = parent.read_arc();
        let edge = edge.read_arc();

        let exp_qt = gtr.expQt(edge.weight().unwrap_or(0.0)).t().to_owned();
        let parent = &parent.dense_partitions[si];
        let edge = &edge.dense_partitions[si];
        let mut dis = Array2::ones((*length, alphabet.n_canonical()));
        let log_lh = parent.msgs_to_children[&name].log_lh;
        if let Some(transmission) = &edge.transmission {
          for r in transmission {
            dis
              .slice_mut(s![r.0..r.1, ..])
              .assign(&(parent.msgs_to_children[&name].dis.slice(s![r.0..r.1, ..]).dot(&exp_qt)));
          }
        } else {
          dis = &dis * &parent.msgs_to_children[&name].dis.dot(&exp_qt);
        }
        msgs_from_parents.insert(name.clone(), DenseSeqDis { dis, log_lh });
      }

      msgs_from_parents.insert(o!("children"), seq_info.msg_to_parents.clone()); // HACK
      seq_info.profile = combine_dense_messages(&msgs_from_parents).unwrap();
      seq_info.seq.sequence = assign_sequence(seq_info, alphabet);
      seq_info.msgs_to_children = btreemap![];
      for cname in seq_info.msgs_from_children.keys() {
        let child_msg = &seq_info.msgs_from_children[cname];
        let mut dis = seq_info.profile.dis.clone();
        dis /= &child_msg.dis;
        let log_lh = seq_info.profile.log_lh - child_msg.log_lh;
        seq_info
          .msgs_to_children
          .insert(cname.clone(), DenseSeqDis { dis, log_lh });
      }
    }

    GraphTraversalContinuation::Continue
  });
}

fn calculate_root_state_dense(graph: &DenseGraph, partitions: &[PartitionLikelihood]) -> f64 {
  let mut total_log_lh = 0.0;
  for root_node in graph.get_roots() {
    let dense_partitions = &mut root_node.write_arc().payload().write_arc().dense_partitions;
    for (si, seq_info) in dense_partitions.iter_mut().enumerate() {
      let PartitionLikelihood { gtr, alphabet, .. } = &partitions[si];

      let mut log_lh = seq_info.msg_to_parents.log_lh;
      let mut dis = &seq_info.msg_to_parents.dis * &gtr.pi;

      let norm = dis.sum_axis(Axis(1));
      log_lh += log(&norm).sum();
      dis = (&dis.t() / &norm).t().to_owned();
      seq_info.profile = DenseSeqDis { dis, log_lh };
      seq_info.seq.sequence = assign_sequence(seq_info, alphabet);
      seq_info.msgs_to_children.clear();
      for cname in seq_info.msgs_from_children.keys() {
        // This division operation can cause 'division by 0' problems. Note that profile = prod(msgs[child], child in children) * msgs_from_parent.
        // Alternatively, msgs_to_children could be calculated as prod(msgs[sibling], sibling in children/child) * msgs_from_parent.
        // This latter redundant calculation is done in the sparse case, but division is more efficient in particular in case of polytomies.
        let child_msg = &seq_info.msgs_from_children[cname];
        let dis = &seq_info.profile.dis / &child_msg.dis;
        let log_lh = seq_info.profile.log_lh - child_msg.log_lh;
        seq_info
          .msgs_to_children
          .insert(cname.clone(), DenseSeqDis { dis, log_lh });
      }

      total_log_lh += seq_info.profile.log_lh;
    }
  }
  total_log_lh
}

pub fn run_marginal_dense(graph: &DenseGraph, partitions: Vec<PartitionLikelihoodWithAln>) -> Result<f64, Report> {
  attach_seqs_to_graph(graph, &partitions)?;

  let partitions = partitions.into_iter().map(PartitionLikelihood::from).collect_vec();

  ingroup_profiles_dense(graph, &partitions);
  let log_lh = calculate_root_state_dense(graph, &partitions);
  outgroup_profiles_dense(graph, &partitions);
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
  use super::*;
  use crate::alphabet::alphabet::AlphabetName;
  use crate::gtr::get_gtr::{jc69, JC69Params};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::json::{json_write_str, JsonPretty};
  use crate::io::nwk::nwk_read_str;
  use crate::utils::string::vec_to_string;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;

  #[test]
  fn test_ancestral_reconstruction_marginal_dense() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(indoc! {r#"
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
    "#})?;

    let expected = read_many_fasta_str(indoc! {r#"
      >root
      TCGGCGCTGTATTGAC
      >AB
      ACATCGCTGTA-TGAC
      >CD
      TCGGCGGTGTATTG--
    "#})?
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
    run_marginal_dense(&graph, partitions)?;

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
}
