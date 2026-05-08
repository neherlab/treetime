#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::optimize_unified::initial_guess_mixed;
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use crate::seq::indel::InDel;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::seq::Seq;

  /// All-zero branch length tree. Auto mode (overwrite_valid=false) treats
  /// zero BL as "valid" and skips the edge - except when indels are present.
  const TREE_ZERO_BL: &str = "((A:0.0,B:0.0)AB:0.0,C:0.0)root:0.0;";

  #[test]
  fn test_initial_guess_auto_preserves_zero_bl_without_indels() -> Result<(), Report> {
    let (graph, partitions) = setup_dense(TREE_ZERO_BL)?;

    initial_guess_mixed(&graph, &partitions, false)?;

    // All edges should remain zero: no indels, zero is valid
    for edge_ref in graph.get_edges() {
      let bl = edge_ref
        .read_arc()
        .payload()
        .read_arc()
        .branch_length()
        .unwrap_or(f64::NAN);
      assert!(bl == 0.0, "Without indels, Auto mode should preserve zero BL, got {bl}");
    }
    Ok(())
  }

  #[test]
  fn test_initial_guess_auto_overrides_zero_bl_with_indels() -> Result<(), Report> {
    let (graph, partitions) = setup_dense(TREE_ZERO_BL)?;

    // Inject an indel on the first edge
    let edge_key = graph.get_edges()[0].read_arc().key();
    {
      let partition = partitions[0].write_arc();
      let edge_data = partition.edges[&edge_key].clone();
      drop(partition);
      let mut partition = partitions[0].write_arc();
      let edge_entry = partition.edges.entry(edge_key).or_insert(edge_data);
      edge_entry.indels.push(InDel {
        range: (4, 7),
        seq: Seq::default(),
        deletion: true,
      });
    }

    initial_guess_mixed(&graph, &partitions, false)?;

    // The indel-bearing edge should now have a positive BL
    let bl = graph.get_edges()[0]
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap_or(0.0);
    assert!(
      bl > 0.0,
      "Auto mode should override zero BL on indel-bearing edge, got {bl}"
    );

    // Non-indel edges should remain zero
    for edge_ref in graph.get_edges().iter().skip(1) {
      let bl = edge_ref
        .read_arc()
        .payload()
        .read_arc()
        .branch_length()
        .unwrap_or(f64::NAN);
      assert!(bl == 0.0, "Non-indel edge should remain zero, got {bl}");
    }
    Ok(())
  }

  mod helpers {
    use super::*;

    pub fn setup_dense(newick: &str) -> Result<(GraphAncestral, Vec<Arc<RwLock<PartitionMarginalDense>>>), Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln = read_many_fasta_str(
        indoc! {r#"
          >A
          AAAACCCCGGGGTTTT
          >B
          CCCCGGGGTTTTAAAA
          >C
          GGGGTTTTAAAACCCC
        "#},
        &alphabet,
      )?;
      let graph: GraphAncestral = nwk_read_str(newick)?;

      let partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense {
        index: 0,
        gtr: jc69(JC69Params::default())?,
        alphabet,
        length: get_common_length(&aln)?,
        nodes: btreemap! {},
        edges: btreemap! {},
      }))];

      initialize_marginal(&graph, &partitions, &aln)?;
      update_marginal(&graph, &partitions)?;

      Ok((graph, partitions))
    }
  }
  use helpers::*;
}
