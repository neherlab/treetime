#[cfg(test)]
mod tests {
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_reconstruction, simple_alignment,
  };
  use crate::optimize::__tests__::test_initial_guess_mode::tests::TREE_ZERO_BL;
  use crate::optimize::__tests__::test_initial_guess_mode::tests::helpers::{
    get_branch_lengths, inject_indel_on_first_edge, setup_dense_with_marginal,
  };
  use crate::optimize::__tests__::test_optimize_indel::tests::{
    inject_indels_on_first_edge, setup_identical_reconstruction,
  };
  use crate::optimize::dispatch::run_optimize_mixed_inner;
  use crate::optimize::gather::{
    gather_edge_contributions, gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts,
  };
  use crate::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
  use crate::optimize::run_loop::{apply_initial_guess_mode, run_optimize_loop};
  use crate::progress::NoopProgress;
  use crate::seq::indel::InDel;
  use crate::test_utils::sparse_reconstruction_mut;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  #[test]
  fn test_no_indels_drops_indel_contribution_from_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_with_names = nwk_parsed.names();
    let mut graph_with = nwk_parsed.graph;
    let mut branch_lengths_with = nwk_parsed.branch_lengths;
    let mut reconstruction_with = setup_reconstruction(&graph_with, &graph_with_names, &aln, &mut branch_lengths_with)?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_without_names = nwk_parsed.names();
    let mut graph_without = nwk_parsed.graph;
    let mut branch_lengths_without = nwk_parsed.branch_lengths;
    let mut reconstruction_without =
      setup_reconstruction(&graph_without, &graph_without_names, &aln, &mut branch_lengths_without)?;

    let first_edge_key = graph_with.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths_with.insert(first_edge_key, Some(0.1));
    sparse_reconstruction_mut(&mut reconstruction_with)
      .partition
      .obs_edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    let first_edge_key_without = graph_without.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths_without.insert(first_edge_key_without, Some(0.1));
    sparse_reconstruction_mut(&mut reconstruction_without)
      .partition
      .obs_edges
      .get_mut(&first_edge_key_without)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    let names_tt_4 = graph_with_names.clone();
    let result_with = run_optimize_loop(
      &mut graph_with,
      reconstruction_with,
      1,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths_with,
      &names_tt_4,
    )?;

    let names_tt_3 = graph_without_names.clone();
    let result_without = run_optimize_loop(
      &mut graph_without,
      reconstruction_without,
      1,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      true,
      TopologyOps::default(),
      branch_lengths_without,
      &names_tt_3,
    )?;

    assert!(
      result_without.lh_history[0] > result_with.lh_history[0],
      "no_indels=true should produce higher likelihood than with indels: {} vs {}",
      result_without.lh_history[0].value(),
      result_with.lh_history[0].value()
    );
    Ok(())
  }

  #[test]
  fn test_no_indels_optimizer_ignores_indel_counts() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let mut reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;

    let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    branch_lengths.insert(first_edge_key, Some(0.05));
    sparse_reconstruction_mut(&mut reconstruction)
      .partition
      .obs_edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    let (reconstruction, _) = reconstruction.marginal_update(&graph, &branch_lengths_or_zero(&branch_lengths))?;

    let bl_before = branch_lengths.clone();

    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
    run_optimize_mixed_inner(
      &graph,
      total_length,
      &contributions,
      &indel_counts,
      BranchOptMethod::BrentSqrt,
      0.0,
      true,
      &mut branch_lengths,
    )?;

    assert_ne!(bl_before, branch_lengths, "Optimizer should modify branch lengths");
    Ok(())
  }

  #[test]
  fn test_no_indels_matches_no_indel_data() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_no_flag_names = nwk_parsed.names();
    let mut graph_no_flag = nwk_parsed.graph;
    let mut branch_lengths_nf = nwk_parsed.branch_lengths;
    let reconstruction_nf = setup_reconstruction(&graph_no_flag, &graph_no_flag_names, &aln, &mut branch_lengths_nf)?;

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_flag_names = nwk_parsed.names();
    let mut graph_flag = nwk_parsed.graph;
    let mut branch_lengths_f = nwk_parsed.branch_lengths;
    let reconstruction_f = setup_reconstruction(&graph_flag, &graph_flag_names, &aln, &mut branch_lengths_f)?;

    let names_tt_2 = graph_no_flag_names.clone();
    let result_no_flag = run_optimize_loop(
      &mut graph_no_flag,
      reconstruction_nf,
      3,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
      TopologyOps::default(),
      branch_lengths_nf,
      &names_tt_2,
    )?;

    let names_tt_1 = graph_flag_names.clone();
    let result_flag = run_optimize_loop(
      &mut graph_flag,
      reconstruction_f,
      3,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      true,
      TopologyOps::default(),
      branch_lengths_f,
      &names_tt_1,
    )?;

    assert_eq!(
      result_no_flag.lh_history.len(),
      result_flag.lh_history.len(),
      "Same number of iterations expected"
    );
    for (i, (lh_nf, lh_f)) in result_no_flag
      .lh_history
      .iter()
      .zip(&result_flag.lh_history)
      .enumerate()
    {
      assert_abs_diff_eq!(lh_nf.value(), lh_f.value(), epsilon = 1e-10);
    }
    Ok(())
  }

  #[test]
  fn test_no_indels_initial_guess_never_accepts_zero_bl_with_indels() -> Result<(), Report> {
    let (graph, names, mut reconstruction, mut branch_lengths) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    inject_indel_on_first_edge(&graph, &mut reconstruction)?;
    let total_length = reconstruction.sequence_length();
    let indel_counts = gather_edge_indel_counts(&graph, &reconstruction);
    let sub_counts = gather_edge_sub_counts(&graph, &reconstruction)?;
    let effective_lengths = gather_edge_effective_lengths(&graph, &reconstruction)?;
    let result = apply_initial_guess_mode(
      &graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      InitialGuessMode::Never,
      true,
      &mut branch_lengths,
      &names,
      &NoopProgress,
    );
    assert!(
      result.is_ok(),
      "no_indels=true should accept zero-BL indel edges in Never mode, got: {result:?}"
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::auto(  InitialGuessMode::Auto)]
  #[case::always(InitialGuessMode::Always)]
  #[trace]
  fn test_no_indels_initial_guess_ignores_indel_counts(#[case] mode: InitialGuessMode,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_with_indel_names = nwk_parsed.names();
    let graph_with_indel = nwk_parsed.graph;
    let mut branch_lengths_with_indel = nwk_parsed.branch_lengths;
    let mut reconstruction_with_indel = setup_identical_reconstruction(&graph_with_indel, &graph_with_indel_names, &mut branch_lengths_with_indel)?;
    let indels = vec![InDel::del((0, 2), Seq::try_from_str("AC")?)?];
    inject_indels_on_first_edge(
      &graph_with_indel,
      &mut reconstruction_with_indel,
      &indels,
    );

    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let graph_without_indel_names = nwk_parsed.names();
    let graph_without_indel = nwk_parsed.graph;
    let mut branch_lengths_without_indel = nwk_parsed.branch_lengths;
    let reconstruction_without_indel = setup_identical_reconstruction(&graph_without_indel, &graph_without_indel_names, &mut branch_lengths_without_indel)?;

    let total_length_with_indel = reconstruction_with_indel.sequence_length();
    let indel_counts_with_indel = gather_edge_indel_counts(&graph_with_indel, &reconstruction_with_indel);
    let sub_counts_with_indel = gather_edge_sub_counts(&graph_with_indel, &reconstruction_with_indel)?;
    let effective_lengths_with_indel = gather_edge_effective_lengths(&graph_with_indel, &reconstruction_with_indel)?;
    apply_initial_guess_mode(&graph_with_indel, total_length_with_indel, &indel_counts_with_indel, &sub_counts_with_indel, &effective_lengths_with_indel, mode, true, &mut branch_lengths_with_indel, &graph_with_indel_names, &NoopProgress)?;
    let total_length_without_indel = reconstruction_without_indel.sequence_length();
    let indel_counts_without_indel = gather_edge_indel_counts(&graph_without_indel, &reconstruction_without_indel);
    let sub_counts_without_indel = gather_edge_sub_counts(&graph_without_indel, &reconstruction_without_indel)?;
    let effective_lengths_without_indel = gather_edge_effective_lengths(&graph_without_indel, &reconstruction_without_indel)?;
    apply_initial_guess_mode(&graph_without_indel, total_length_without_indel, &indel_counts_without_indel, &sub_counts_without_indel, &effective_lengths_without_indel, mode, true, &mut branch_lengths_without_indel, &graph_without_indel_names, &NoopProgress)?;

    let expected = get_branch_lengths(&graph_without_indel, &branch_lengths_without_indel);
    let actual = get_branch_lengths(&graph_with_indel, &branch_lengths_with_indel);
    assert_eq!(expected, actual);
    Ok(())
  }
}
