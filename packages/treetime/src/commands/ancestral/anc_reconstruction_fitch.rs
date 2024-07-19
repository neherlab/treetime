use crate::commands::ancestral::anc_graph::{Edge, Node};
use crate::graph::graph::{Graph, SafeNode};
use crate::seq::representation::{apply_non_nuc_changes_inplace, pre_order_intrusive};
use eyre::Report;
use std::borrow::Cow;

/// Reconstruct ancestral sequences using Fitch parsimony.
///
/// Calls visitor function for every ancestral node, providing the node itself and its reconstructed sequence.
/// Optionally reconstructs leaf sequences.
pub fn ancestral_reconstruction_fitch(
  graph: &Graph<Node, Edge>,
  include_leaves: bool,
  mut visitor: impl FnMut(&Node, &[char]),
) -> Result<(), Report> {
  let root = graph.get_exactly_one_root()?;
  let mut root_seq = { root.read_arc().payload().read_arc().seq.clone() };

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
      let seq = if include_leaves {
        let mut seq = seq.to_owned();
        apply_non_nuc_changes_inplace(&node, &mut seq);
        Cow::from(seq)
      } else {
        Cow::from(seq)
      };

      visitor(&node, &seq);
    },
  );

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::ancestral::anc_reconstruction_fitch::ancestral_reconstruction_fitch;
  use crate::io::json::{json_write_str, JsonPretty};
  use crate::io::nwk::nwk_read_str;
  use crate::o;
  use crate::seq::representation::compress_sequences;
  use crate::utils::random::get_random_number_generator;
  use crate::utils::string::vec_to_string;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  #[rstest]
  fn test_seq_ancestral_reconstruction_with_leaves() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let mut rng = get_random_number_generator(Some(42));

    let inputs = BTreeMap::from([
      (o!("A"), o!("ACATCGCCNNA--G")),
      (o!("B"), o!("GCATCCCTGTA-NG")),
      (o!("C"), o!("CCGGCGATGTATTG")),
      (o!("D"), o!("TCGGCCGTGTRTTG")),
    ]);

    #[rustfmt::skip]
    let expected = BTreeMap::from([
      (o!("A"),    o!("ACATCGCCNNA--G")),
      (o!("AB"),   o!("GCATCGCTGTATTG")),
      (o!("B"),    o!("GCATCCCTGTA-NG")),
      (o!("C"),    o!("CCGGCGATGTATTG")),
      (o!("CD"),   o!("TCGGCGGTGTATTG")),
      (o!("D"),    o!("TCGGCCGTGTRTTG")),
      (o!("root"), o!("TCGGCGGTGTATTG")),
    ]);

    let graph = nwk_read_str::<Node, Edge>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    compress_sequences(&inputs, &graph, &mut rng).unwrap();

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_fitch(&graph, true, |node, seq| {
      actual.insert(node.name.clone(), vec_to_string(seq.to_owned()));
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }
}
