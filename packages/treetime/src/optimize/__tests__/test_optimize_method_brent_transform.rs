#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {

  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_partitions, simple_alignment,
  };
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_indel_counts, total_sequence_length};

  use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
  use crate::optimize::method_brent::{brent_bracket, brent_log_inner, brent_sqrt_inner};

  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::zero_boundary::min_branch_length_for_indels;

  use eyre::Report;

  use rstest::rstest;

  use treetime_graph::graph::Graph;

  use treetime_io::nwk::nwk_read_str;

  use crate::optimize::__tests__::test_optimize_method::tests::helpers::*;

  #[test]
  fn test_optimize_method_brent_sqrt_transform_round_trip() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln = simple_alignment()?;
    let (dense_mixed_partitions, sparse_mixed_partitions) =
      setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;
    let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    let mut contributions_by_edge =
      gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let contributions = contributions_by_edge
      .remove(&edge_key)
      .expect("first edge present in gathered contributions");

    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let one_mutation = 1.0 / total_length as f64;
    let branch_length = 0.01;

    let result = brent_sqrt_inner(branch_length, &contributions, 0, 0.0, 0.0, one_mutation).unwrap();
    assert!(
      result >= 0.0,
      "brent_sqrt_inner result must be non-negative, got {result}"
    );
    assert!(
      result.is_finite(),
      "brent_sqrt_inner result must be finite, got {result}"
    );

    let lh_opt = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result)
      .expect("valid branch length")
      .value();
    if result > 1e-10 {
      let lh_below = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result * 0.99)
        .expect("valid branch length")
        .value();
      let lh_above = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result * 1.01)
        .expect("valid branch length")
        .value();
      assert!(
        lh_opt >= lh_below - 1e-10,
        "sqrt: lh at opt ({lh_opt}) < lh below ({lh_below})"
      );
      assert!(
        lh_opt >= lh_above - 1e-10,
        "sqrt: lh at opt ({lh_opt}) < lh above ({lh_above})"
      );
    }

    Ok(())
  }

  #[test]
  fn test_optimize_method_brent_log_transform_round_trip() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let aln = simple_alignment()?;
    let (dense_mixed_partitions, sparse_mixed_partitions) =
      setup_partitions(&graph, &names, &aln, &mut branch_lengths)?;
    let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    let mut contributions_by_edge =
      gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let contributions = contributions_by_edge
      .remove(&edge_key)
      .expect("first edge present in gathered contributions");

    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let one_mutation = 1.0 / total_length as f64;
    let branch_length = 0.01;

    let result = brent_log_inner(branch_length, &contributions, 0, 0.0, 0.0, one_mutation).unwrap();
    assert!(result > 0.0, "brent_log_inner result must be positive, got {result}");
    assert!(
      result.is_finite(),
      "brent_log_inner result must be finite, got {result}"
    );

    let lh_opt = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result)
      .expect("valid branch length")
      .value();
    let lh_below = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result * 0.99)
      .expect("valid branch length")
      .value();
    let lh_above = evaluate_with_indels_log_lh_only(&contributions, 0, 0.0, result * 1.01)
      .expect("valid branch length")
      .value();
    assert!(
      lh_opt >= lh_below - 1e-10,
      "log: lh at opt ({lh_opt}) < lh below ({lh_below})"
    );
    assert!(
      lh_opt >= lh_above - 1e-10,
      "log: lh at opt ({lh_opt}) < lh above ({lh_above})"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::brent(     BranchOptMethod::Brent)]
  #[case::brent_sqrt(BranchOptMethod::BrentSqrt)]
  #[case::brent_log( BranchOptMethod::BrentLog)]
  #[trace]
  fn test_optimize_method_brent_bracket_validity(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_mixed_partitions, sparse_mixed_partitions, indel_rate) = setup_with_indels(&graph, &names, &mut branch_lengths, 4)?;
    let total_length = total_sequence_length(&dense_mixed_partitions, &sparse_mixed_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_mixed_partitions, &sparse_mixed_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_mixed_partitions, &sparse_mixed_partitions);

    let input_bl = branch_lengths[&graph.get_edges().collect::<Vec<_>>()[0].key()].unwrap_or(0.0);
    let one_mutation = 1.0 / total_length as f64;
    let min_bl = min_branch_length_for_indels(4, one_mutation);
    let (lower, upper) = brent_bracket(input_bl, min_bl, one_mutation);

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;

    let bl = first_edge_bl(&graph, &branch_lengths);
    let lh_opt = eval_combined_first_edge(&graph, &dense_mixed_partitions, &sparse_mixed_partitions, indel_rate, bl)?;

    let lh_lower = eval_combined_first_edge(&graph, &dense_mixed_partitions, &sparse_mixed_partitions, indel_rate, lower)?;
    let lh_upper = eval_combined_first_edge(&graph, &dense_mixed_partitions, &sparse_mixed_partitions, indel_rate, upper)?;

    assert!(
      lh_opt >= lh_lower - 1e-10,
      "{method:?} optimum lh ({lh_opt}) < lower bracket lh ({lh_lower})"
    );
    assert!(
      lh_opt >= lh_upper - 1e-10,
      "{method:?} optimum lh ({lh_opt}) < upper bracket lh ({lh_upper})"
    );

    Ok(())
  }
}
