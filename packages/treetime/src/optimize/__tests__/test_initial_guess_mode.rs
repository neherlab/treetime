#[cfg(test)]
pub(super) mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{JC69Params, jc69};
  use crate::optimize::branch_length::invalid_branch_length_descriptions;
  use crate::optimize::dispatch::guess_initial_branch_lengths;
  use crate::optimize::gather::InitialGuessEvidence;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::params::ExistingBranchLengths;
  use crate::optimize::params::InitialGuessMode;
  use crate::optimize::run_loop::{apply_initial_guess_mode, invalid_branch_length_warning};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;
  use crate::partition::marginal::reconstruction::DenseReconstruction;
  use crate::partition::set::PartitionSet;
  use crate::progress::NoopProgress;
  use crate::test_utils::leaf_seq_inputs;
  use crate::test_utils::{deletion, dense_reconstruction_mut};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node_names::NodeNames;

  use treetime_io::fasta::{FastaRecord, fasta_read};
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;
  use treetime_primitives::Seq;

  const TREE_WITH_LENGTHS: &str = "((A:0.1,B:0.2)AB:0.05,C:0.3)root:0.01;";

  const TREE_WITHOUT_LENGTHS: &str = "((A,B)AB,C)root;";

  pub(crate) const TREE_ZERO_BL: &str = "((A:0.0,B:0.0)AB:0.0,C:0.0)root:0.0;";

  #[test]
  fn test_initial_guess_mode_default_is_auto() {
    assert_eq!(InitialGuessMode::Auto, InitialGuessMode::default());
  }

  #[test]
  fn test_initial_guess_mode_detects_nan_from_newick() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_WITHOUT_LENGTHS.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    assert!(!invalid_branch_length_descriptions(&graph, &branch_lengths, &names)?.is_empty());
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_detects_explicit_nan() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_WITH_LENGTHS.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let nan_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths.insert(nan_edge_key, Some(f64::NAN));
    assert!(!invalid_branch_length_descriptions(&graph, &branch_lengths, &names)?.is_empty());
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_detects_negative_branch_length() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_WITH_LENGTHS.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    set_first_branch_length(&graph, &mut branch_lengths, -0.1);
    assert!(!invalid_branch_length_descriptions(&graph, &branch_lengths, &names)?.is_empty());
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_describes_all_invalid_edges_and_warning() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_WITH_LENGTHS.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
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
    let nwk_parsed = nwk_read(TREE_WITH_LENGTHS.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    assert!(invalid_branch_length_descriptions(&graph, &branch_lengths, &names)?.is_empty());
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_preserves_valid_lengths() -> Result<(), Report> {
    let (graph, names, reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let before = get_branch_lengths(&graph, &branch_lengths);

    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    guess_initial_branch_lengths(
      &graph,
      &indels,
      &evidence,
      ExistingBranchLengths::Keep,
      &mut branch_lengths,
    )?;

    let after = get_branch_lengths(&graph, &branch_lengths);
    assert_eq!(
      before, after,
      "Auto mode must not overwrite valid finite branch lengths"
    );
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_fills_none_from_newick() -> Result<(), Report> {
    let (graph, names, reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_WITHOUT_LENGTHS)?;

    for edge_ref in graph.get_edges() {
      let bl = branch_lengths[&edge_ref.key()];
      assert!(
        bl.is_none(),
        "Expected None from newick parser for missing branch lengths"
      );
    }

    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    guess_initial_branch_lengths(
      &graph,
      &indels,
      &evidence,
      ExistingBranchLengths::Keep,
      &mut branch_lengths,
    )?;

    for edge_ref in graph.get_edges() {
      let bl = branch_lengths[&edge_ref.key()];
      assert!(
        bl.is_some_and(|v| v.is_finite() && v >= 0.0),
        "Expected finite non-negative branch length, got {bl:?}"
      );
    }
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_auto_fills_only_missing_edges() -> Result<(), Report> {
    let (graph, names, reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let original_lengths = get_branch_lengths(&graph, &branch_lengths);

    let nan_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths.insert(nan_edge_key, Some(f64::NAN));

    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    guess_initial_branch_lengths(
      &graph,
      &indels,
      &evidence,
      ExistingBranchLengths::Keep,
      &mut branch_lengths,
    )?;

    let updated_lengths = get_branch_lengths(&graph, &branch_lengths);

    assert!(
      updated_lengths[0].is_finite() && updated_lengths[0] >= 0.0,
      "Missing edge should be filled with finite value"
    );

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
    let (graph, names, mut reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    set_first_branch_length(&graph, &mut branch_lengths, -0.1);
    if has_indels {
      inject_indel_on_first_edge(&graph, &mut reconstruction)?;
    }

    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    let result = apply_initial_guess_mode(
      &graph,
      &indels,
      &evidence,
      mode,
      &mut branch_lengths,
      &names,
      &NoopProgress,
    );

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
    let (graph, names, reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let before = get_branch_lengths(&graph, &branch_lengths);

    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    guess_initial_branch_lengths(
      &graph,
      &indels,
      &evidence,
      ExistingBranchLengths::Overwrite,
      &mut branch_lengths,
    )?;

    let after = get_branch_lengths(&graph, &branch_lengths);
    assert_ne!(before, after, "Always mode must overwrite existing branch lengths");
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_accepts_complete_tree() -> Result<(), Report> {
    let (graph, names, reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    let before = get_branch_lengths(&graph, &branch_lengths);
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    apply_initial_guess_mode(
      &graph,
      &indels,
      &evidence,
      InitialGuessMode::Never,
      &mut branch_lengths,
      &names,
      &NoopProgress,
    )?;
    let after = get_branch_lengths(&graph, &branch_lengths);
    assert_eq!(before, after, "Never mode must leave branch lengths unchanged");
    Ok(())
  }

  #[test]
  fn test_initial_guess_mode_never_rejects_nan_tree() -> Result<(), Report> {
    let (graph, names, reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_WITHOUT_LENGTHS)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    let result = apply_initial_guess_mode(
      &graph,
      &indels,
      &evidence,
      InitialGuessMode::Never,
      &mut branch_lengths,
      &names,
      &NoopProgress,
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
    let (graph, names, reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    let before = get_branch_lengths(&graph, &branch_lengths);
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    apply_initial_guess_mode(
      &graph,
      &indels,
      &evidence,
      InitialGuessMode::Never,
      &mut branch_lengths,
      &names,
      &NoopProgress,
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
    let (graph, names, mut reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    inject_indel_on_first_edge(&graph, &mut reconstruction)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    let result = apply_initial_guess_mode(
      &graph,
      &indels,
      &evidence,
      InitialGuessMode::Never,
      &mut branch_lengths,
      &names,
      &NoopProgress,
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
    let (graph, names, mut reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_WITH_LENGTHS)?;
    inject_indel_on_first_edge(&graph, &mut reconstruction)?;
    let before = get_branch_lengths(&graph, &branch_lengths);
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
    let evidence = InitialGuessEvidence::gather(&graph, &reconstruction)?;
    apply_initial_guess_mode(
      &graph,
      &indels,
      &evidence,
      InitialGuessMode::Never,
      &mut branch_lengths,
      &names,
      &NoopProgress,
    )?;
    let after = get_branch_lengths(&graph, &branch_lengths);
    assert_eq!(
      before, after,
      "Never mode must accept positive branch lengths on indel-bearing edges"
    );
    Ok(())
  }

  pub(crate) mod helpers {
    use super::*;

    pub(crate) fn get_branch_lengths(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> Vec<f64> {
      graph
        .get_edges()
        .map(|e| branch_lengths[&e.key()].unwrap_or(f64::NAN))
        .collect()
    }

    pub(crate) fn set_first_branch_length(
      graph: &Graph,
      branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
      branch_length: f64,
    ) {
      let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
      branch_lengths.insert(edge_key, Some(branch_length));
    }

    pub(crate) fn set_branch_length_by_target_name(
      names: &NodeNames,
      graph: &Graph,
      branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
      target_name: &str,
      branch_length: f64,
    ) {
      let edge = graph
        .get_edges()
        .into_iter()
        .find(|edge_ref| {
          let target = edge_ref.target();
          graph
            .get_node(target)
            .is_some_and(|node| &names[node.key()] == target_name)
        })
        .expect("target node must have an incoming edge");
      let edge_key = edge.key();
      branch_lengths.insert(edge_key, Some(branch_length));
    }

    pub(crate) fn inject_indel_on_first_edge(graph: &Graph, reconstruction: &mut PartitionSet) -> Result<(), Report> {
      let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
      dense_reconstruction_mut(reconstruction)
        .edges
        .estimates
        .get_mut(&edge_key)
        .unwrap()
        .indels = vec![deletion((4, 7), Seq::try_from_str("ACG")?)];
      Ok(())
    }

    pub(crate) fn setup_dense_with_marginal(
      newick: &str,
    ) -> Result<(Graph, NodeNames, PartitionSet, BTreeMap<GraphEdgeKey, Option<f64>>), Report> {
      let alphabet = Alphabet::new(AlphabetName::Nuc)?;
      let aln: Vec<AlignmentRecord> = test_alignment()?.into_iter().map(AlignmentRecord::from).collect();
      let nwk_parsed = nwk_read(newick.as_bytes())?;
      let names = nwk_parsed.names.clone();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;

      let partition = PartitionMarginalDense::new(alphabet, &graph, &leaf_seq_inputs(&graph, &names, aln))?;
      let reconstruction = PartitionSet::Dense(vec![DenseReconstruction::seeded(
        partition,
        jc69(JC69Params::default())?,
      )]);

      let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

      Ok((graph, names, reconstruction, branch_lengths))
    }

    fn test_alignment() -> Result<Vec<FastaRecord>, Report> {
      let alphabet = Alphabet::default();
      fasta_read(
        indoc! {r#"
          >A
          AAAACCCCGGGGTTTT
          >B
          CCCCGGGGTTTTAAAA
          >C
          GGGGTTTTAAAACCCC
        "#}
        .as_bytes(),
        &alphabet,
      )
    }
  }
  use helpers::*;
}
