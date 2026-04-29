#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::commands::ancestral::fitch::get_common_length;
  use crate::commands::ancestral::marginal::{initialize_marginal, update_marginal};
  use crate::commands::optimize::args::InitialGuessMode;
  use crate::commands::optimize::optimize_unified::initial_guess_mixed;
  use crate::commands::optimize::run::{
    any_edge_missing_branch_length, any_indel_edge_has_zero_branch_length, apply_initial_guess_mode,
  };
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::representation::partition::marginal_dense::PartitionMarginalDense;
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  const TREE_WITH_LENGTHS: &str = "((A:0.1,B:0.2)AB:0.05,C:0.3)root:0.01;";

  const TREE_WITHOUT_LENGTHS: &str = "((A,B)AB,C)root;";

  const TREE_ZERO_BL: &str = "((A:0.0,B:0.0)AB:0.0,C:0.0)root:0.0;";

  #[test]
  fn test_initial_guess_mode_default_is_auto() {
    assert_eq!(InitialGuessMode::Auto, InitialGuessMode::default());
  }

  #[test]
  fn test_initial_guess_mode_detects_nan_from_newick() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITHOUT_LENGTHS)?;
    assert!(any_edge_missing_branch_length(&graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_detects_explicit_nan() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITH_LENGTHS)?;
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(f64::NAN));
    assert!(any_edge_missing_branch_length(&graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_no_missing_when_all_finite() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_WITH_LENGTHS)?;
    assert!(!any_edge_missing_branch_length(&graph));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_preserves_valid_lengths() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let before = get_branch_lengths(&graph);

    initial_guess_mixed(&graph, &partitions, false)?;

    let after = get_branch_lengths(&graph);
    assert_eq!(
      before, after,
      "Auto mode must not overwrite valid finite branch lengths"
    );
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_fills_nan_from_newick() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_WITHOUT_LENGTHS)?;

    // Before: all edges have NaN (from bio crate newick parser)
    for edge_ref in graph.get_edges() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length();
      assert!(bl.is_some_and(|v| v.is_nan()), "Expected NaN from newick parser");
    }

    initial_guess_mixed(&graph, &partitions, false)?;

    // After: all edges have finite non-negative values
    for edge_ref in graph.get_edges() {
      let bl = edge_ref.read_arc().payload().read_arc().branch_length();
      assert!(
        bl.is_some_and(|v| v.is_finite() && v >= 0.0),
        "Expected finite non-negative branch length, got {bl:?}"
      );
    }
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_fills_only_missing_edges() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let original_lengths = get_branch_lengths(&graph);

    // Set one edge to NaN (simulating a missing branch length)
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(f64::NAN));

    initial_guess_mixed(&graph, &partitions, false)?;

    let updated_lengths = get_branch_lengths(&graph);

    // The previously-NaN edge now has a finite value
    assert!(
      updated_lengths[0].is_finite() && updated_lengths[0] >= 0.0,
      "Missing edge should be filled with finite value"
    );

    // All other edges retain their original values
    for i in 1..original_lengths.len() {
      assert_abs_diff_eq!(original_lengths[i], updated_lengths[i], epsilon = 1e-15);
    }
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_always_overwrites_all() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let before = get_branch_lengths(&graph);

    initial_guess_mixed(&graph, &partitions, true)?;

    let after = get_branch_lengths(&graph);
    assert_ne!(before, after, "Always mode must overwrite existing branch lengths");
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_accepts_complete_tree() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let before = get_branch_lengths(&graph);
    apply_initial_guess_mode(&graph, &partitions, InitialGuessMode::Never, false)?;
    let after = get_branch_lengths(&graph);
    assert_eq!(before, after, "Never mode must leave branch lengths unchanged");
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_rejects_nan_tree() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_WITHOUT_LENGTHS)?;
    let result = apply_initial_guess_mode(&graph, &partitions, InitialGuessMode::Never, false);
    let err = result.expect_err("Never mode must reject a tree with NaN branch lengths");
    let msg = format!("{err:?}");
    assert!(
      msg.contains("--branch-length-initial-guess=never"),
      "Never-mode error must mention the offending flag, got: {msg}"
    );
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_accepts_zero_bl_without_indels() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    let before = get_branch_lengths(&graph);
    apply_initial_guess_mode(&graph, &partitions, InitialGuessMode::Never, false)?;
    let after = get_branch_lengths(&graph);
    assert_eq!(
      before, after,
      "Never mode must accept all-zero branch lengths when no indels are present"
    );
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_rejects_zero_bl_with_indels() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    let result = apply_initial_guess_mode(&graph, &partitions, InitialGuessMode::Never, false);
    let err = result.expect_err("Never mode must reject zero branch length on an indel-bearing edge");
    let msg = format!("{err:?}");
    assert!(
      msg.contains("--branch-length-initial-guess=never"),
      "Never-mode error must mention the offending flag, got: {msg}"
    );
    assert!(
      msg.contains("indel"),
      "Never-mode error must explain the indel-specific rejection, got: {msg}"
    );
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_accepts_positive_bl_with_indels() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    let before = get_branch_lengths(&graph);
    apply_initial_guess_mode(&graph, &partitions, InitialGuessMode::Never, false)?;
    let after = get_branch_lengths(&graph);
    assert_eq!(
      before, after,
      "Never mode must accept positive branch lengths on indel-bearing edges"
    );
    Ok(())
  }

  #[test]
  fn test_any_indel_edge_has_zero_bl_false_without_indels() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    assert!(
      !any_indel_edge_has_zero_branch_length(&graph, &partitions),
      "Without any indels, no indel-bearing zero-BL edges should be detected"
    );
    Ok(())
  }

  #[test]
  fn test_any_indel_edge_has_zero_bl_false_with_positive_bl() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    assert!(
      !any_indel_edge_has_zero_branch_length(&graph, &partitions),
      "Positive branch length on indel-bearing edge must not trigger the zero-BL check"
    );
    Ok(())
  }

  #[test]
  fn test_any_indel_edge_has_zero_bl_true_with_indel_and_zero_bl() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    assert!(
      any_indel_edge_has_zero_branch_length(&graph, &partitions),
      "Zero branch length on an indel-bearing edge must be detected"
    );
    Ok(())
  }

  mod helpers {
    use super::*;

    pub fn get_branch_lengths(graph: &GraphAncestral) -> Vec<f64> {
      graph
        .get_edges()
        .iter()
        .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(f64::NAN))
        .collect()
    }

    /// Attach a single 3-base deletion indel to the partition's entry for
    /// the first graph edge. Marginal initialization populates the edge
    /// map, so the entry always exists by the time this helper is called.
    pub fn inject_indel_on_first_edge(
      graph: &GraphAncestral,
      partitions: &[Arc<RwLock<PartitionMarginalDense>>],
    ) -> Result<(), Report> {
      let edge_key = graph.get_edges()[0].read_arc().key();
      for partition in partitions {
        let mut partition = partition.write_arc();
        partition.edges.get_mut(&edge_key).unwrap().indels = vec![InDel::del((4, 7), Seq::try_from_str("ACG")?)];
      }
      Ok(())
    }

    pub fn setup_dense_with_marginal(
      newick: &str,
    ) -> Result<(GraphAncestral, Vec<Arc<RwLock<PartitionMarginalDense>>>), Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln = test_alignment()?;
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

    fn test_alignment() -> Result<Vec<FastaRecord>, Report> {
      let alphabet = Alphabet::default();
      read_many_fasta_str(
        indoc! {r#"
          >A
          AAAACCCCGGGGTTTT
          >B
          CCCCGGGGTTTTAAAA
          >C
          GGGGTTTTAAAACCCC
        "#},
        &alphabet,
      )
    }
  }
  use helpers::*;
}
