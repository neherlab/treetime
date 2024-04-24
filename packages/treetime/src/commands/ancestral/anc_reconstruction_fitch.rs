use crate::commands::ancestral::anc_graph::{Edge, Node};
use crate::graph::breadth_first::{GraphTraversalContinuation, GraphTraversalPhase};
use crate::graph::graph::{Graph, GraphNodeForward};
use crate::seq::representation::{apply_muts_inplace, apply_non_nuc_changes_inplace};
use eyre::Report;

/// Reconstruct ancestral sequences using Fitch parsimony.
///
/// Calls visitor function for every ancestral node, providing the node itself and its reconstructed sequence.
/// Optionally reconstructs leaf sequences.
pub fn ancestral_reconstruction_fitch(
  graph: &Graph<Node, Edge>,
  include_leaves: bool,
  visitor: impl Fn(&Node, &[char]) + Send + Sync,
) -> Result<(), Report> {
  graph.par_iter_breadth_first_biphasic_forward(
    |(mut node, phase): (GraphNodeForward<Node, Edge>, GraphTraversalPhase)| {
      match phase {
        GraphTraversalPhase::Pre => {
          if !include_leaves && node.is_leaf {
            return GraphTraversalContinuation::Continue;
          }

          // Apply node mutations to the parent sequence, and remember the result
          if let Some((parent, _)) = node.get_one_parent() {
            let mut seq = { parent.read_arc().seq.clone() };
            apply_muts_inplace(&node.payload, &mut seq);
            if include_leaves && node.is_leaf {
              apply_non_nuc_changes_inplace(&node.payload, &mut seq);
            }
            node.payload.seq = seq;
          }

          visitor(&node.payload, &node.payload.seq);
        }
        GraphTraversalPhase::Post => {
          // Second visit, after all children are visited. Deallocate sequences - they are no longer needed.
          if !node.is_root {
            node.payload.seq = vec![];
          }
        }
      }
      GraphTraversalContinuation::Continue
    },
  );

  Ok(())
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::commands::ancestral::anc_reconstruction_fitch::ancestral_reconstruction_fitch;
  use crate::graph::create_graph_from_nwk::create_graph_from_nwk_str;
  use crate::io::json::{json_stringify, JsonPretty};
  use crate::o;
  use crate::seq::representation::compress_sequences;
  use crate::utils::random::get_random_number_generator;
  use crate::utils::string::vec_to_string;
  use eyre::Report;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;

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

    let graph = create_graph_from_nwk_str::<Node, Edge>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    compress_sequences(&inputs, &graph, &mut rng).unwrap();

    let actual = Arc::new(RwLock::new(BTreeMap::new()));
    ancestral_reconstruction_fitch(&graph, true, |node, seq| {
      actual
        .write_arc()
        .insert(node.name.clone(), vec_to_string(seq.to_owned()));
    })?;

    assert_eq!(
      json_stringify(&expected, JsonPretty(false))?,
      json_stringify(&actual, JsonPretty(false))?
    );

    Ok(())
  }
}
