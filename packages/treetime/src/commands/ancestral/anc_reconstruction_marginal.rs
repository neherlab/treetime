use crate::commands::ancestral::anc_graph::{Edge, Node};
use crate::commands::ancestral::outgroup_profiles::outgroup_profiles;
use crate::commands::ancestral::subtree_profiles::subtree_profiles;
use crate::graph::graph::{Graph, SafeNode};
use crate::gtr::gtr::GTR;
use crate::seq::representation::{apply_non_nuc_changes_inplace, post_order_intrusive, pre_order_intrusive};
use eyre::Report;
use itertools::Itertools;
use ndarray_stats::QuantileExt;

/// Reconstruct ancestral sequences using marginal method.
///
/// Calls visitor function for every ancestral node, providing the node itself and its reconstructed sequence.
/// Optionally reconstructs leaf sequences.
pub fn ancestral_reconstruction_marginal(
  graph: &Graph<Node, Edge>,
  gtr: &GTR,
  include_leaves: bool,
  mut visitor: impl FnMut(&Node, &[char]),
) -> Result<(), Report> {
  let root = graph.get_exactly_one_root()?;
  let mut root_seq = { root.read_arc().payload().read_arc().seq.clone() };
  let mut logLH = 0.0;

  post_order_intrusive(graph, &root, None, &mut root_seq, &mut |node, edge, seq| {
    let node = node.write_arc();

    let children: Vec<SafeNode<Node>> = graph
      .children_of(&node)
      .into_iter()
      .map(|(child, _)| child)
      .collect_vec();

    let is_leaf = node.is_leaf();
    let mut node = node.payload().write_arc();
    let edge = edge.map(|e| e.payload().read_arc());
    subtree_profiles(
      graph,
      &mut node,
      is_leaf,
      edge.as_deref(),
      &children,
      seq,
      gtr,
      &mut logLH,
    );
  });

  pre_order_intrusive(graph, &root, &mut root_seq, &mut |node, seq| {
    let node = node.write_arc();

    let parent = graph
      .one_parent_of(&node)
      .unwrap()
      .map(|parent| parent.read_arc().payload().read_arc());

    let mut node = node.payload().write_arc();
    outgroup_profiles(&mut node, parent.as_deref(), seq, &mut logLH, gtr);
  });

  pre_order_intrusive(
    graph,
    &root,
    &mut root_seq,
    &mut |node: &SafeNode<Node>, seq: &[char]| {
      let node = node.read_arc();

      if !include_leaves && node.is_leaf() {
        return;
      }

      let node = node.payload().read_arc();
      let mut seq = seq.to_owned();
      for (&pos, vec) in &node.profile_variable {
        seq[pos] = gtr.alphabet.char(vec.argmax().unwrap());
      }

      if include_leaves {
        apply_non_nuc_changes_inplace(&node, &mut seq);
      }

      visitor(&node, &seq);
    },
  );

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::anc_graph::{Edge, Node};
  use crate::gtr::gtr::GTRParams;
  use crate::io::json::{json_write_str, JsonPretty};
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use crate::seq::representation::compress_sequences;
  use crate::utils::random::get_random_number_generator;
  use crate::utils::string::vec_to_string;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  #[rstest]
  fn test_ancestral_reconstruction_marginal() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let mut rng = get_random_number_generator(Some(42));

    let inputs = BTreeMap::from([
      (o!("A"), o!("ACATCGCCNNA--G")),
      (o!("B"), o!("GCATCCCTGTA-NG")),
      (o!("C"), o!("CCGGCGATGTATTG")),
      (o!("D"), o!("TCGGCCGTGTRTTG")),
    ]);

    let graph = nwk_read_str::<Node, Edge>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    compress_sequences(&inputs, &graph, &mut rng).unwrap();

    let alphabet = Alphabet::new(AlphabetName::NucNogap)?;

    let gtr = GTR::new(GTRParams {
      alphabet,
      mu: 1.0,
      W: None,
      pi: array![0.2, 0.3, 0.15, 0.45],
    })?;

    let mut actual = btreemap! {};
    ancestral_reconstruction_marginal(&graph, &gtr, true, |node, seq| {
      actual.insert(node.name.clone(), vec_to_string(seq.to_owned()));
    })?;

    #[rustfmt::skip]
    let expected = BTreeMap::from([
      (o!("A"),    o!("ACATCGCCNNA--G")),
      (o!("AB"),   o!("ACATCGCTGTATTG")),
      (o!("B"),    o!("GCATCCCTGTA-NG")),
      (o!("C"),    o!("CCGGCGATGTATTG")),
      (o!("CD"),   o!("TCGGCGGTGTATTG")),
      (o!("D"),    o!("TCGGCCGTGTRTTG")),
      (o!("root"), o!("TCGGCGCTGTATTG")),
    ]);

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }
}
