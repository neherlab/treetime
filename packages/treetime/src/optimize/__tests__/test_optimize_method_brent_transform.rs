#[cfg(test)]
mod tests {

  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_reconstruction, simple_alignment,
  };
  use crate::optimize::__tests__::test_optimize_method::tests::helpers::*;
  use crate::optimize::dispatch::run_optimize_mixed;
  use crate::optimize::gather::gather_edge_contributions;
  use crate::optimize::indel::{IndelPolicy, IndelTerm};
  use crate::optimize::likelihood::evaluate_with_indels_log_lh_only;
  use crate::optimize::method_brent::{BranchLengthBracket, brent_bracket, brent_sqrt_inner};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::zero_boundary::min_branch_length_for_indels;
  use eyre::Report;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::cast::ToF64Exact;

  #[test]
  fn test_optimize_method_brent_sqrt_transform_round_trip() -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let aln = simple_alignment()?;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;
    let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    let mut contributions_by_edge = gather_edge_contributions(&graph, &reconstruction)?;
    let contributions = contributions_by_edge
      .remove(&edge_key)
      .expect("first edge present in gathered contributions");

    let total_length = reconstruction.sequence_length();
    let one_mutation = 1.0 / total_length.to_f64_exact();
    let branch_length = 0.01;

    let result = brent_sqrt_inner(branch_length, &contributions, IndelTerm::new(0, 0.0), 0.0, one_mutation).unwrap();
    assert!(
      result >= 0.0,
      "brent_sqrt_inner result must be non-negative, got {result}"
    );
    assert!(
      result.is_finite(),
      "brent_sqrt_inner result must be finite, got {result}"
    );

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
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let (reconstruction, indel_rate) = setup_with_indels(&graph, &names, &mut branch_lengths, 4)?;
    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);

    let input_bl = branch_lengths[&graph.get_edges().collect::<Vec<_>>()[0].key()].unwrap_or(0.0);
    let one_mutation = 1.0 / total_length.to_f64_exact();
    let min_bl = min_branch_length_for_indels(IndelTerm::new(4, indel_rate), one_mutation);
    let BranchLengthBracket { lower, upper } = brent_bracket(input_bl, min_bl, one_mutation);

    run_optimize_mixed(
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
}
