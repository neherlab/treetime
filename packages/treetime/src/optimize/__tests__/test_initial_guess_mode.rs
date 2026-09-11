#[cfg(test)]
pub mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::marginal::{initialize_marginal, marginal_update, profile_branch_lengths};
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::branch_length::invalid_branch_length_descriptions;
  use crate::optimize::dispatch::initial_guess_mixed;
  use crate::optimize::params::InitialGuessMode;
  use crate::optimize::run_loop::{
    any_indel_edge_has_zero_branch_length, apply_initial_guess_mode, invalid_branch_length_warning,
  };
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::alignment::get_common_length;
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use treetime_graph::node::GraphNodeKey;

  use parking_lot::RwLock;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::sync::Arc;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::{FastaRecord, read_many_fasta_str};
  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_primitives::Seq;

  const TREE_WITH_LENGTHS: &str = "((A:0.1,B:0.2)AB:0.05,C:0.3)root:0.01;";

  pub const TREE_WITHOUT_LENGTHS: &str = "((A,B)AB,C)root;";

  pub const TREE_ZERO_BL: &str = "((A:0.0,B:0.0)AB:0.0,C:0.0)root:0.0;";

  #[test]
  fn test_initial_guess_mode_default_is_auto() {
    assert_eq!(InitialGuessMode::Auto, InitialGuessMode::default());
  }

  #[test]
  fn test_initial_guess_mode_detects_nan_from_newick() -> Result<(), Report> {
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(TREE_WITHOUT_LENGTHS)?;
    let graph: GraphAncestral = graph;
    assert!(!invalid_branch_length_descriptions(&graph, &branch_lengths, &names)?.is_empty());
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_detects_explicit_nan() -> Result<(), Report> {
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_WITH_LENGTHS)?;
    let graph: GraphAncestral = graph;
    let nan_edge_key = graph.get_edges()[0].read_arc().key();
    branch_lengths.insert(nan_edge_key, Some(f64::NAN));
    assert!(!invalid_branch_length_descriptions(&graph, &branch_lengths, &names)?.is_empty());
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_detects_negative_branch_length() -> Result<(), Report> {
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_WITH_LENGTHS)?;
    let graph: GraphAncestral = graph;
    set_first_branch_length(&graph, &mut branch_lengths, -0.1);
    assert!(!invalid_branch_length_descriptions(&graph, &branch_lengths, &names)?.is_empty());
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_describes_all_invalid_edges_and_warning() -> Result<(), Report> {
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_WITH_LENGTHS)?;
    let graph: GraphAncestral = graph;
    set_branch_length_by_target_name(&names, &graph, &mut branch_lengths, "A", -0.1);
    set_branch_length_by_target_name(&names, &graph, &mut branch_lengths, "C", f64::INFINITY);

    let descriptions = invalid_branch_length_descriptions(&graph, &branch_lengths, &names)?;
    let expected = vec!["AB -> A: -0.1", "root -> C: inf"];
    assert_eq!(expected, descriptions);

    let warning = invalid_branch_length_warning(&descriptions).expect("invalid edges produce a warning");
    assert!(warning.contains("AB -> A: -0.1"));
    assert!(warning.contains("root -> C: inf"));
    assert!(warning.contains("--branch-length-initial-guess=always"));
    assert_eq!(None, invalid_branch_length_warning(&[]));
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_no_missing_when_all_finite() -> Result<(), Report> {
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str(TREE_WITH_LENGTHS)?;
    let graph: GraphAncestral = graph;
    assert!(invalid_branch_length_descriptions(&graph, &branch_lengths, &names)?.is_empty());
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_preserves_valid_lengths() -> Result<(), Report> {
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let before = get_branch_lengths(&graph, &branch_lengths);

    initial_guess_mixed(&graph, &partitions, false, false, &mut branch_lengths)?;

    let after = get_branch_lengths(&graph, &branch_lengths);
    assert_eq!(
      before, after,
      "Auto mode must not overwrite valid finite branch lengths"
    );
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_fills_none_from_newick() -> Result<(), Report> {
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_WITHOUT_LENGTHS)?;

    // Before: all edges have None (pest parser returns None for missing branch lengths)
    for edge_ref in graph.get_edges() {
      let bl = branch_lengths[&edge_ref.read_arc().key()];
      assert!(
        bl.is_none(),
        "Expected None from newick parser for missing branch lengths"
      );
    }

    initial_guess_mixed(&graph, &partitions, false, false, &mut branch_lengths)?;

    // After: all edges have finite non-negative values
    for edge_ref in graph.get_edges() {
      let bl = branch_lengths[&edge_ref.read_arc().key()];
      assert!(
        bl.is_some_and(|v| v.is_finite() && v >= 0.0),
        "Expected finite non-negative branch length, got {bl:?}"
      );
    }
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_fills_only_missing_edges() -> Result<(), Report> {
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let original_lengths = get_branch_lengths(&graph, &branch_lengths);

    // Set one edge to NaN (simulating a missing branch length)
    let nan_edge_key = graph.get_edges()[0].read_arc().key();
    branch_lengths.insert(nan_edge_key, Some(f64::NAN));

    initial_guess_mixed(&graph, &partitions, false, false, &mut branch_lengths)?;

    let updated_lengths = get_branch_lengths(&graph, &branch_lengths);

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

  #[rustfmt::skip]
  #[rstest]
  #[case::auto_without_indels(  InitialGuessMode::Auto,   false, false)]
  #[case::auto_with_indels(     InitialGuessMode::Auto,   true,  false)]
  #[case::always_without_indels(InitialGuessMode::Always, false, false)]
  #[case::always_with_indels(   InitialGuessMode::Always, true,  false)]
  #[case::never_without_indels( InitialGuessMode::Never,  false, true)]
  #[case::never_with_indels(    InitialGuessMode::Never,  true,  true)]
  #[trace]
  fn test_initial_guess_mode_negative_branch_length_matrix(
    #[case] mode: InitialGuessMode,
    #[case] has_indels: bool,
    #[case] expects_error: bool,
  ) -> Result<(), Report> {
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    set_first_branch_length(&graph, &mut branch_lengths, -0.1);
    if has_indels {
      inject_indel_on_first_edge(&graph, &partitions)?;
    }

    let result = apply_initial_guess_mode(&graph, &partitions, mode, false, &mut branch_lengths, &names);

    if expects_error {
      let error = result.expect_err("Never mode must reject a negative branch length");
      let message = format!("{error:?}");
      assert!(message.contains("-0.1"), "Error must include the invalid value, got: {message}");
      return Ok(());
    }
    result?;

    let branch_length = get_branch_lengths(&graph, &branch_lengths)[0];
    assert!(
      branch_length.is_finite() && branch_length >= 0.0,
      "Initial-guess mode must replace a negative branch length, got {branch_length}"
    );
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_always_overwrites_all() -> Result<(), Report> {
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let before = get_branch_lengths(&graph, &branch_lengths);

    initial_guess_mixed(&graph, &partitions, true, false, &mut branch_lengths)?;

    let after = get_branch_lengths(&graph, &branch_lengths);
    assert_ne!(before, after, "Always mode must overwrite existing branch lengths");
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_accepts_complete_tree() -> Result<(), Report> {
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let before = get_branch_lengths(&graph, &branch_lengths);
    apply_initial_guess_mode(
      &graph,
      &partitions,
      InitialGuessMode::Never,
      false,
      &mut branch_lengths,
      &names,
    )?;
    let after = get_branch_lengths(&graph, &branch_lengths);
    assert_eq!(before, after, "Never mode must leave branch lengths unchanged");
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_rejects_nan_tree() -> Result<(), Report> {
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_WITHOUT_LENGTHS)?;
    let result = apply_initial_guess_mode(
      &graph,
      &partitions,
      InitialGuessMode::Never,
      false,
      &mut branch_lengths,
      &names,
    );
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
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    let before = get_branch_lengths(&graph, &branch_lengths);
    apply_initial_guess_mode(
      &graph,
      &partitions,
      InitialGuessMode::Never,
      false,
      &mut branch_lengths,
      &names,
    )?;
    let after = get_branch_lengths(&graph, &branch_lengths);
    assert_eq!(
      before, after,
      "Never mode must accept all-zero branch lengths when no indels are present"
    );
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_rejects_zero_bl_with_indels() -> Result<(), Report> {
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    let result = apply_initial_guess_mode(
      &graph,
      &partitions,
      InitialGuessMode::Never,
      false,
      &mut branch_lengths,
      &names,
    );
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
    let (graph, names, partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    let before = get_branch_lengths(&graph, &branch_lengths);
    apply_initial_guess_mode(
      &graph,
      &partitions,
      InitialGuessMode::Never,
      false,
      &mut branch_lengths,
      &names,
    )?;
    let after = get_branch_lengths(&graph, &branch_lengths);
    assert_eq!(
      before, after,
      "Never mode must accept positive branch lengths on indel-bearing edges"
    );
    Ok(())
  }

  #[test]
  fn test_any_indel_edge_has_zero_bl_false_without_indels() -> Result<(), Report> {
    let (graph, names, partitions, branch_lengths) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    assert!(
      !any_indel_edge_has_zero_branch_length(&graph, &partitions, &branch_lengths),
      "Without any indels, no indel-bearing zero-BL edges should be detected"
    );
    Ok(())
  }

  #[test]
  fn test_any_indel_edge_has_zero_bl_false_with_positive_bl() -> Result<(), Report> {
    let (graph, names, partitions, branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    assert!(
      !any_indel_edge_has_zero_branch_length(&graph, &partitions, &branch_lengths),
      "Positive branch length on indel-bearing edge must not trigger the zero-BL check"
    );
    Ok(())
  }

  #[test]
  fn test_any_indel_edge_has_zero_bl_true_with_indel_and_zero_bl() -> Result<(), Report> {
    let (graph, names, partitions, branch_lengths) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    assert!(
      any_indel_edge_has_zero_branch_length(&graph, &partitions, &branch_lengths),
      "Zero branch length on an indel-bearing edge must be detected"
    );
    Ok(())
  }

  pub mod helpers {
    use super::*;

    pub fn get_branch_lengths(
      graph: &GraphAncestral,
      branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
    ) -> Vec<f64> {
      graph
        .get_edges()
        .iter()
        .map(|e| branch_lengths[&e.read_arc().key()].unwrap_or(f64::NAN))
        .collect()
    }

    pub fn set_first_branch_length(
      graph: &GraphAncestral,
      branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
      branch_length: f64,
    ) {
      let edge_key = graph.get_edges()[0].read_arc().key();
      branch_lengths.insert(edge_key, Some(branch_length));
    }

    pub fn set_branch_length_by_target_name(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &GraphAncestral,
      branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
      target_name: &str,
      branch_length: f64,
    ) {
      let edge = graph
        .get_edges()
        .into_iter()
        .find(|edge_ref| {
          let target = edge_ref.read_arc().target();
          graph
            .get_node(target)
            .is_some_and(|node| names.get(&node.read_arc().key()).and_then(|x| x.as_deref()) == Some(target_name))
        })
        .expect("target node must have an incoming edge");
      let edge_key = edge.read_arc().key();
      branch_lengths.insert(edge_key, Some(branch_length));
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
        partition.data.edges.get_mut(&edge_key).unwrap().indels = vec![InDel::del((4, 7), Seq::try_from_str("ACG")?)?];
      }
      Ok(())
    }

    pub fn setup_dense_with_marginal(
      newick: &str,
    ) -> Result<
      (
        GraphAncestral,
        BTreeMap<GraphNodeKey, Option<String>>,
        Vec<Arc<RwLock<PartitionMarginalDense>>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln = test_alignment()?;
      let NwkParse {
        graph,
        names,
        branch_lengths,
        ..
      } = nwk_read_str(newick)?;
      let graph: GraphAncestral = graph;

      let partitions = vec![Arc::new(RwLock::new(PartitionMarginalDense::new(
        0,
        jc69(JC69Params::default())?,
        alphabet,
        get_common_length(&aln)?,
      )))];

      initialize_marginal(
        &graph,
        &profile_branch_lengths(&branch_lengths),
        &partitions,
        &aln,
        &names,
      )?
      .value();
      marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &partitions)?.value();

      Ok((graph, names, partitions, branch_lengths))
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
