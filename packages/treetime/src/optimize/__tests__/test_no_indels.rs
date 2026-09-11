#[cfg(test)]
mod tests {
  use crate::ancestral::marginal::{marginal_update, profile_branch_lengths};
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
  use crate::optimize::params::{BranchOptMethod, InitialGuessMode, TopologyOps};
  use crate::optimize::run_loop::apply_initial_guess_mode;
  use crate::optimize::run_loop::optimize_partition_view;
  use crate::optimize::run_loop::run_optimize_loop;
  use crate::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use rstest::rstest;
  use treetime_io::nwk::{NwkParse, nwk_read_str};
  use treetime_primitives::Seq;

  #[test]
  fn test_no_indels_drops_indel_contribution_from_likelihood() -> Result<(), Report> {
    let aln = simple_alignment()?;
    let NwkParse {
      graph: mut graph_with,
      names: graph_with_names,
      branch_lengths: mut branch_lengths_with,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let (mut dense_with, mut sparse_with) =
      setup_partitions(&graph_with, &graph_with_names, &aln, &mut branch_lengths_with)?;
    let mixed_with = optimize_partition_view(&dense_with, &sparse_with);

    let NwkParse {
      graph: mut graph_without,
      names: graph_without_names,
      branch_lengths: mut branch_lengths_without,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let (mut dense_without, mut sparse_without) =
      setup_partitions(&graph_without, &graph_without_names, &aln, &mut branch_lengths_without)?;
    let mixed_without = optimize_partition_view(&dense_without, &sparse_without);

    let first_edge_key = graph_with.get_edges()[0].read_arc().key();
    branch_lengths_with.insert(first_edge_key, Some(0.1));
    sparse_with[0].edges.get_mut(&first_edge_key).unwrap().indels =
      vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    let first_edge_key_without = graph_without.get_edges()[0].read_arc().key();
    branch_lengths_without.insert(first_edge_key_without, Some(0.1));
    sparse_without[0].edges.get_mut(&first_edge_key_without).unwrap().indels =
      vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    let names_tt_4 = graph_with_names.clone();
    let result_with = run_optimize_loop(
      &mut graph_with,
      &mut sparse_with,
      &mut dense_with,
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
      &mut sparse_without,
      &mut dense_without,
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
    let NwkParse {
      graph,
      names,
      mut branch_lengths,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let graph: GraphAncestral = graph;
    let (dense_partitions, mut sparse_partitions) = setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;

    let first_edge_key = graph.get_edges()[0].read_arc().key();
    branch_lengths.insert(first_edge_key, Some(0.05));
    sparse_partitions[0].edges.get_mut(&first_edge_key).unwrap().indels =
      vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)?];

    marginal_update(&graph, &profile_branch_lengths(&branch_lengths), &mut sparse_partitions)?.value();

    let bl_before = branch_lengths.clone();

    run_optimize_mixed_inner(
      &graph,
      &optimize_partition_view(&dense_partitions, &sparse_partitions),
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
    let NwkParse {
      graph: mut graph_no_flag,
      names: graph_no_flag_names,
      branch_lengths: mut branch_lengths_nf,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let (mut dense_nf, mut sparse_nf) =
      setup_partitions(&graph_no_flag, &graph_no_flag_names, &aln, &mut branch_lengths_nf)?;
    let mixed_nf = optimize_partition_view(&dense_nf, &sparse_nf);

    let NwkParse {
      graph: mut graph_flag,
      names: graph_flag_names,
      branch_lengths: mut branch_lengths_f,
      ..
    } = nwk_read_str(TREE_NEWICK)?;
    let (mut dense_f, mut sparse_f) = setup_partitions(&graph_flag, &graph_flag_names, &aln, &mut branch_lengths_f)?;
    let mixed_f = optimize_partition_view(&dense_f, &sparse_f);

    let names_tt_2 = graph_no_flag_names.clone();
    let result_no_flag = run_optimize_loop(
      &mut graph_no_flag,
      &mut sparse_nf,
      &mut dense_nf,
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
      &mut sparse_f,
      &mut dense_f,
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
    let (graph, names, mut partitions, mut branch_lengths) = setup_dense_with_marginal(TREE_ZERO_BL)?;
    inject_indel_on_first_edge(&graph, &mut partitions)?;
    let result = apply_initial_guess_mode(
      &graph,
      &optimize_partition_view(&partitions, &[]),
      InitialGuessMode::Never,
      true,
      &mut branch_lengths,
      &names,
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
    let NwkParse { graph: graph_with_indel, names: graph_with_indel_names, branch_lengths: mut branch_lengths_with_indel, .. } = nwk_read_str(TREE_NEWICK)?;
    let (mut dense_with_indel, mut sparse_with_indel) = setup_identical_partitions(&graph_with_indel, &graph_with_indel_names, &mut branch_lengths_with_indel)?;
    let indels = vec![InDel::del((0, 2), Seq::try_from_str("AC")?)?];
    inject_indels_on_first_edge(
      &graph_with_indel,
      &mut dense_with_indel,      &mut sparse_with_indel,      &indels,
    );

    let NwkParse { graph: graph_without_indel, names: graph_without_indel_names, branch_lengths: mut branch_lengths_without_indel, .. } = nwk_read_str(TREE_NEWICK)?;
    let (mut dense_without_indel, mut sparse_without_indel) = setup_identical_partitions(&graph_without_indel, &graph_without_indel_names, &mut branch_lengths_without_indel)?;

    apply_initial_guess_mode(&graph_with_indel, &optimize_partition_view(&dense_with_indel, &sparse_with_indel), mode, true, &mut branch_lengths_with_indel, &graph_with_indel_names)?;
    apply_initial_guess_mode(&graph_without_indel, &optimize_partition_view(&dense_without_indel, &sparse_without_indel), mode, true, &mut branch_lengths_without_indel, &graph_without_indel_names)?;

    let expected = get_branch_lengths(&graph_without_indel, &branch_lengths_without_indel);
    let actual = get_branch_lengths(&graph_with_indel, &branch_lengths_with_indel);
    assert_eq!(expected, actual);
    Ok(())
  }
}
