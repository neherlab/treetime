#[cfg(test)]
mod tests {

  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_reconstruction, simple_alignment,
  };
  use crate::optimize::__tests__::test_optimize_method::tests::helpers::*;
  use crate::optimize::dispatch::optimize_edge_lengths;
  use crate::optimize::gather::gather_edge_contributions;
  use crate::optimize::indel::{IndelPolicy, IndelTerm};
  use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
  use crate::optimize::method_brent::{BranchLengthBracket, brent_bracket, brent_sqrt_minimize, solve_brent};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::solver_outcome::Convergence;
  use crate::optimize::zero_boundary::min_branch_length_for_indels;
  use eyre::Report;
  use helpers::QuadraticInSqrt;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::cast::ToF64Exact;

  #[test]
  fn test_optimize_method_brent_sqrt_transform_round_trip() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let aln = simple_alignment()?;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;
    let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    let mut contributions_by_edge = gather_edge_contributions(&graph, &reconstruction)?;
    let contributions = contributions_by_edge
      .remove(&edge_key)
      .expect("first edge present in gathered contributions");

    let total_length = reconstruction.total_length();
    let one_mutation = 1.0 / total_length.to_f64_exact();
    let branch_length = 0.01;

    let solution = solve_brent(branch_length, &contributions, IndelTerm::new(0, 0.0), 0.0, one_mutation).unwrap();
    assert_eq!(Convergence::Converged, solution.convergence);
    let result = solution.branch_length;
    assert!(result >= 0.0, "solve_brent result must be non-negative, got {result}");
    assert!(result.is_finite(), "solve_brent result must be finite, got {result}");

    let lh_opt = evaluate_with_indels_log_lh_only(&contributions, IndelTerm::new(0, 0.0), result)
      .expect("valid branch length")
      .value();
    if result > 1e-10 {
      let lh_below = evaluate_with_indels_log_lh_only(&contributions, IndelTerm::new(0, 0.0), result * 0.99)
        .expect("valid branch length")
        .value();
      let lh_above = evaluate_with_indels_log_lh_only(&contributions, IndelTerm::new(0, 0.0), result * 1.01)
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
  fn test_optimize_method_brent_sqrt_bracket_validity() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let (reconstruction, indel_rate) = setup_with_indels(&graph, &names, &mut branch_lengths, 4)?;
    let total_length = reconstruction.total_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths)?;

    let input_bl = branch_lengths[&graph.get_edges().collect::<Vec<_>>()[0].key()].unwrap_or(0.0);
    let one_mutation = 1.0 / total_length.to_f64_exact();
    let min_bl = min_branch_length_for_indels(IndelTerm::new(4, indel_rate), one_mutation);
    let BranchLengthBracket { lower, upper } = brent_bracket(input_bl, min_bl, one_mutation);

    optimize_edge_lengths(
      &graph,
      total_length,
      &contributions,
      &indels,
      BranchOptMethod::Brent,
      &mut branch_lengths,
    )?;

    let bl = first_edge_bl(&graph, &branch_lengths);
    let lh_opt = eval_combined_first_edge(&graph, &reconstruction, indel_rate, bl)?;

    let lh_lower = eval_combined_first_edge(&graph, &reconstruction, indel_rate, lower)?;
    let lh_upper = eval_combined_first_edge(&graph, &reconstruction, indel_rate, upper)?;

    assert!(
      lh_opt >= lh_lower - 1e-10,
      "Brent optimum lh ({lh_opt}) < lower bracket lh ({lh_lower})"
    );
    assert!(
      lh_opt >= lh_upper - 1e-10,
      "Brent optimum lh ({lh_opt}) < upper bracket lh ({lh_upper})"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_iteration(  1, Convergence::IterationLimit)]
  #[case::default_budget(50, Convergence::Converged)]
  #[trace]
  fn test_optimize_method_brent_sqrt_minimize_reports_termination_reason(
    #[case] max_iters: u64,
    #[case] expected: Convergence,
  ) -> Result<(), Report> {
    let bracket = BranchLengthBracket { lower: 1e-4, upper: 1.0 };
    let solution = brent_sqrt_minimize(&QuadraticInSqrt { minimum: 0.3 }, &bracket, max_iters)?;
    assert_eq!(expected, solution.convergence);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::minimum_inside_bracket(0.3,   Convergence::Converged  )]
  #[case::minimum_below_bracket( 0.001, Convergence::Converged  )]
  #[case::minimum_above_bracket( 2.0,   Convergence::SearchBound)]
  #[trace]
  fn test_optimize_method_brent_sqrt_minimize_reports_search_bound_only_at_upper_end(
    #[case] minimum: f64,
    #[case] expected: Convergence,
  ) -> Result<(), Report> {
    let bracket = BranchLengthBracket { lower: 1e-4, upper: 1.0 };
    let solution = brent_sqrt_minimize(&QuadraticInSqrt { minimum }, &bracket, 50)?;
    assert_eq!(
      expected, solution.convergence,
      "the bracket is s in [0.01, 1]; only a minimum of (s - {minimum})^2 above s = 1 leaves the cost still falling at the upper end"
    );
    Ok(())
  }

  mod helpers {
    use argmin::core::{CostFunction, Error};

    pub(super) struct QuadraticInSqrt {
      pub(super) minimum: f64,
    }

    impl CostFunction for &QuadraticInSqrt {
      type Param = f64;
      type Output = f64;

      fn cost(&self, s: &Self::Param) -> Result<Self::Output, Error> {
        Ok((s - self.minimum).powi(2))
      }
    }
  }
}
