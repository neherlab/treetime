#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::update_marginal;
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::optimize::__tests__::test_initial_guess_mode::tests::TREE_ZERO_BL;
  use crate::optimize::__tests__::test_initial_guess_mode::tests::helpers::{
    get_branch_lengths, inject_indel_on_first_edge, setup_dense_with_marginal,
  };
  use crate::optimize::__tests__::test_optimize_indel::tests::{
    inject_indels_on_first_edge, setup_identical_partitions,
  };
  use crate::optimize::dispatch::run_optimize_mixed_inner;
  use crate::optimize::params::{BranchOptMethod, InitialGuessMode};
  use crate::optimize::run_loop::apply_initial_guess_mode;
  use crate::optimize::run_loop::run_optimize_loop;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use rstest::rstest;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  #[test]
  fn test_no_indels_drops_indel_contribution_from_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph_with: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_with, sparse_with, mixed_with) = setup_partitions(&graph_with, &aln)?;

    let mut graph_without: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_without, sparse_without, mixed_without) = setup_partitions(&graph_without, &aln)?;

    let first_edge_key = graph_with.get_edges()[0].read_arc().key();
    graph_with.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(0.1));
    sparse_with[0]
      .write_arc()
      .edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    let first_edge_key_without = graph_without.get_edges()[0].read_arc().key();
    graph_without.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(0.1));
    sparse_without[0]
      .write_arc()
      .edges
      .get_mut(&first_edge_key_without)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    let result_with = run_optimize_loop(
      &mut graph_with,
      &sparse_with,
      &dense_with,
      &mixed_with,
      1,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
    )?;

    let result_without = run_optimize_loop(
      &mut graph_without,
      &sparse_without,
      &dense_without,
      &mixed_without,
      1,
      0.0,
      0.75,
      BranchOptMethod::BrentSqrt,
      true,
    )?;

    assert!(
      result_without.lh_history[0] > result_with.lh_history[0],
      "no_indels=true should produce higher likelihood than with indels: {} vs {}",
      result_without.lh_history[0],
      result_with.lh_history[0]
    );
    Ok(())
  }

  #[test]
  fn test_no_indels_optimizer_ignores_indel_counts() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (_, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let first_edge_key = graph.get_edges()[0].read_arc().key();
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(0.05));
    sparse_partitions[0]
      .write_arc()
      .edges
      .get_mut(&first_edge_key)
      .unwrap()
      .indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    update_marginal(&graph, &sparse_partitions)?;

    let bl_before: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .collect();

    run_optimize_mixed_inner(&graph, &mixed_partitions, BranchOptMethod::BrentSqrt, 0.0, true)?;

    let bl_after: Vec<f64> = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .collect();

    assert_ne!(bl_before, bl_after, "Optimizer should modify branch lengths");
    Ok(())
  }

  #[test]
  fn test_no_indels_matches_no_indel_data() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let mut graph_no_flag: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_nf, sparse_nf, mixed_nf) = setup_partitions(&graph_no_flag, &aln)?;

    let mut graph_flag: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_f, sparse_f, mixed_f) = setup_partitions(&graph_flag, &aln)?;

    let result_no_flag = run_optimize_loop(
      &mut graph_no_flag,
      &sparse_nf,
      &dense_nf,
      &mixed_nf,
      3,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      false,
    )?;

    let result_flag = run_optimize_loop(
      &mut graph_flag,
      &sparse_f,
      &dense_f,
      &mixed_f,
      3,
      0.1,
      0.75,
      BranchOptMethod::BrentSqrt,
      true,
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
      assert_abs_diff_eq!(lh_nf, lh_f, epsilon = 1e-10);
    }
    Ok(())
  }

  #[test]
  fn test_no_indels_initial_guess_never_accepts_zero_bl_with_indels() -> Result<(), Report> {
    let (graph, partitions) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    inject_indel_on_first_edge(&graph, &partitions)?;
    let result = apply_initial_guess_mode(&graph, &partitions, InitialGuessMode::Never, true);
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
  fn test_no_indels_initial_guess_ignores_indel_counts(
    #[case] mode: InitialGuessMode,
  ) -> Result<(), Report> {
    let graph_with_indel: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (dense_with_indel, sparse_with_indel, partitions_with_indel) =
      setup_identical_partitions(&graph_with_indel)?;
    let indels = vec![InDel::del((0, 2), Seq::try_from_str("AC")?)?];
    inject_indels_on_first_edge(
      &graph_with_indel,
      &dense_with_indel,
      &sparse_with_indel,
      &indels,
    );

    let graph_without_indel: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let (_, _, partitions_without_indel) = setup_identical_partitions(&graph_without_indel)?;

    apply_initial_guess_mode(&graph_with_indel, &partitions_with_indel, mode, true)?;
    apply_initial_guess_mode(&graph_without_indel, &partitions_without_indel, mode, true)?;

    let expected = get_branch_lengths(&graph_without_indel);
    let actual = get_branch_lengths(&graph_with_indel);
    assert_eq!(expected, actual);
    Ok(())
  }
}
