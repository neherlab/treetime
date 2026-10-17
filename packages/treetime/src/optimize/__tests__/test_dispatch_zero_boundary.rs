#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::branch_lengths::branch_lengths_or_zero;
  use crate::gtr::get_gtr::{GtrModelName, JC69Params, K80Params, get_gtr_by_name, jc69, k80};
  use crate::gtr::gtr::{GTR, Mu};
  use crate::optimize::dispatch::{guess_initial_branch_lengths, optimize_edge_lengths};
  use crate::optimize::gather::{gather_edge_contributions, gather_edge_effective_lengths, gather_edge_sub_counts};
  use crate::optimize::indel::{IndelPolicy, IndelTerm};
  use crate::optimize::likelihood::evaluate_partitions_log_lh_only;
  use crate::optimize::method_newton::solve_newton;
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::params::ExistingBranchLengths;
  use crate::optimize::run_loop::find_zero_optimal_internal_edges;
  use crate::optimize::solver_outcome::Convergence;
  use crate::optimize::zero_boundary::{ZeroLength, is_zero_branch_optimal, reconcile_zero_boundary};
  use crate::partition::fitch::passes::create_fitch_partition;
  use crate::partition::marginal::reconstruction::{MarginalReconstruction, SparseReconstruction};
  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use crate::test_utils::leaf_seq_inputs;
  use crate::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};
  use eyre::Report;
  use indoc::indoc;
  use ndarray::{Array2, array};
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::fasta::fasta_read;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::AlignmentRecord;

  const IDENTICAL_TREE_NEWICK: &str = "((A:0.1,B:0.1)AB:0.1,(C:0.1,D:0.1)CD:0.1)root:0.1;";

  const IDENTICAL_ALIGNMENT: &str = indoc! {r#"
    >A
    ACGTACGTACGTACGT
    >B
    ACGTACGTACGTACGT
    >C
    ACGTACGTACGTACGT
    >D
    ACGTACGTACGTACGT
  "#};

  fn setup_identical_reconstruction(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    model: GtrModelName,
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<MarginalReconstruction, Report> {
    let aln: Vec<AlignmentRecord> = fasta_read(IDENTICAL_ALIGNMENT.as_bytes(), &Alphabet::default())?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let fitch = create_fitch_partition(
      graph,
      Alphabet::new(AlphabetName::Nuc)?,
      leaf_seq_inputs(graph, names, aln),
    )?;
    let (sparse_partition, sparse_node_states) = fitch.into_marginal_sparse(graph)?;
    let reconstruction = MarginalReconstruction::Sparse(SparseReconstruction::seeded(
      sparse_partition,
      get_gtr_by_name(model)?,
      sparse_node_states,
    ));
    let (reconstruction, _) = reconstruction.marginal_update(graph, &branch_lengths_or_zero(branch_lengths))?;

    let total_length = reconstruction.sequence_length();
    let indels = IndelPolicy::Include.edge_indels(graph, &reconstruction, branch_lengths);
    let sub_counts = gather_edge_sub_counts(graph, &reconstruction)?;
    let effective_lengths = gather_edge_effective_lengths(graph, &reconstruction)?;
    guess_initial_branch_lengths(
      graph,
      total_length,
      &indels,
      &sub_counts,
      &effective_lengths,
      ExistingBranchLengths::Keep,
      branch_lengths,
    )?;

    Ok(reconstruction)
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(BranchOptMethod::Newton)]
  #[case::brent( BranchOptMethod::Brent)]
  #[trace]
  fn test_dispatch_zero_boundary_k80_identical_sequences(#[case] method: BranchOptMethod,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(IDENTICAL_TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_identical_reconstruction(&graph, &names, GtrModelName::K80, &mut branch_lengths)?;
    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);

    for (i, edge_ref) in graph.get_edges().enumerate() {
      let bl = branch_lengths[&edge_ref.key()].unwrap();
      assert!(
        bl > 0.0,
        "precondition: edge {i} must start with positive BL, got {bl}"
      );
    }

    optimize_edge_lengths(&graph, total_length, &contributions, &indels, method, &mut branch_lengths)?;

    for (i, edge_ref) in graph.get_edges().enumerate() {
      let bl = branch_lengths[&edge_ref.key()].unwrap();
      assert!(
        bl == 0.0,
        "{method:?}: edge {i} branch length must be exactly 0 after optimization, got {bl}"
      );
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k80(  GtrModelName::K80)]
  #[case::hky85(GtrModelName::HKY85)]
  #[case::t92(  GtrModelName::T92)]
  #[trace]
  fn test_dispatch_zero_boundary_uncertified_models_all_reach_zero(#[case] model: GtrModelName,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(IDENTICAL_TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_identical_reconstruction(&graph, &names, model, &mut branch_lengths)?;
    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);

    assert_eq!(
      None,
      reconstruction.gtr().certified_decay_rate(),
      "precondition: {model:?} must have no certified single decay rate"
    );

    optimize_edge_lengths(
      &graph,
      total_length,
      &contributions,
      &indels,
      BranchOptMethod::Brent,
      &mut branch_lengths,
    )?;

    for (i, edge_ref) in graph.get_edges().enumerate() {
      let bl = branch_lengths[&edge_ref.key()].unwrap();
      assert!(
        bl == 0.0,
        "{model:?}: edge {i} branch length must be exactly 0 after optimization, got {bl}"
      );
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::jc69(GtrModelName::JC69)]
  #[case::tn93(GtrModelName::TN93)]
  #[trace]
  fn test_dispatch_zero_boundary_certified_models_pre_dispatch_shortcut_reaches_zero(
    #[case] model: GtrModelName,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(IDENTICAL_TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_identical_reconstruction(&graph, &names, model, &mut branch_lengths)?;
    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);

    assert!(
      reconstruction.gtr().certified_decay_rate().is_some(),
      "precondition: {model:?} must have a certified single decay rate"
    );

    let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    let first_edge_contributions = &contributions[&first_edge_key];
    let zero_length = ZeroLength::of_edge(first_edge_contributions, IndelTerm::new(0, 0.0));
    assert!(
      is_zero_branch_optimal(first_edge_contributions, zero_length),
      "precondition: {model:?} identical-sequence contributions must trigger the pre-dispatch shortcut"
    );

    optimize_edge_lengths(
      &graph,
      total_length,
      &contributions,
      &indels,
      BranchOptMethod::Brent,
      &mut branch_lengths,
    )?;
    for (i, edge_ref) in graph.get_edges().enumerate() {
      let bl = branch_lengths[&edge_ref.key()].unwrap();
      assert!(
        bl == 0.0,
        "{model:?}: edge {i} branch length must be exactly 0 after optimization, got {bl}"
      );
    }

    Ok(())
  }

  fn make_dinh_matsen_k80_contribution() -> OptimizationContribution {
    #[rustfmt::skip]
    let coefficients = array![
      [-0.04545042,  0.02261158,  0.02261158, 0.25],
      [ 0.04456328, -0.02228164, -0.02228164, 0.25],
    ];
    OptimizationContribution::Dense(optimize::dense::PartitionContribution {
      coefficients,
      rate_eigvals: array![-1.0, -0.5, -0.5, 0.0],
      gtr: k80(K80Params::default()).unwrap(),
    })
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_positive_candidate_worse_than_zero_returns_grid_point_beating_zero() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let one_mutation = 0.01;
    let start = 0.2;
    let candidate = 1.0;
    let lh_zero = evaluate_partitions_log_lh_only(&contributions, 0.0).unwrap().value();
    let lh_candidate = evaluate_partitions_log_lh_only(&contributions, candidate)
      .unwrap()
      .value();
    assert!(
      lh_candidate < lh_zero,
      "precondition: candidate log_lh({candidate})={lh_candidate} must be worse than log_lh(0)={lh_zero}"
    );

    let result = reconcile(candidate, start, &contributions, IndelTerm::new(0, 0.0), one_mutation);

    let lh_result = evaluate_partitions_log_lh_only(&contributions, result).unwrap().value();
    assert!(
      lh_result > lh_zero,
      "the grid point must beat zero: log_lh({result})={lh_result} vs log_lh(0)={lh_zero}"
    );
    pretty_assert_ulps_eq!(
      0.16152813746142794,
      result,
      max_ulps = 4,
      "expected the best of 100 log-spaced points on [0.001, 0.5] (index 81), next to the local maximum t = 0.16657 of the exact fixture coefficients (numpy evaluation); the global maximum t = 6.37 lies outside the grid"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_degenerate_site_positive_candidate_passes_through() {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.0, 0.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr));
    let contributions = [contribution];

    assert!(
      !contributions[0].all_sites_valid_at_zero(),
      "precondition: degenerate contribution must fail all_sites_valid_at_zero"
    );

    let candidate = 0.01;
    let result = reconcile(candidate, 0.1, &contributions, IndelTerm::new(0, 0.0), 0.001);
    assert!(
      result.to_bits() == candidate.to_bits(),
      "reconcile_zero_boundary must return positive candidate unchanged when a site is degenerate, got {result}"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_degenerate_site_zero_candidate_routes_to_grid() {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.0, 0.0, 0.0, 0.0], [0.5, 0.3, 0.1, 0.1]];
    let contribution = OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr));
    let contributions = [contribution];

    assert!(
      !contributions[0].all_sites_valid_at_zero(),
      "precondition: degenerate contribution must fail all_sites_valid_at_zero"
    );

    let candidate = 0.0;
    let one_mutation = 0.001;
    let result = reconcile(candidate, 0.1, &contributions, IndelTerm::new(0, 0.0), one_mutation);
    assert!(
      result > 0.0,
      "reconcile_zero_boundary must NOT preserve 0.0 when a site is degenerate at zero, got {result}"
    );
    assert!(result.is_finite(), "reconciled result must be finite, got {result}");
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_indel_count_positive_passes_through() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let candidate = 0.005;
    let indel_count = 1;
    let indel_rate = 44.4;
    let one_mutation = 0.01;

    let result = reconcile(
      candidate,
      0.2,
      &contributions,
      IndelTerm::new(indel_count, indel_rate),
      one_mutation,
    );
    assert!(
      result.to_bits() == candidate.to_bits(),
      "reconcile_zero_boundary must return candidate unchanged when indel_count > 0, got {result}"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::brent(BranchOptMethod::Brent)]
  #[trace]
  fn test_dispatch_zero_boundary_topology_cleanup_collects_k80_internal_edges(#[case] method: BranchOptMethod,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read(IDENTICAL_TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let reconstruction = setup_identical_reconstruction(&graph, &names, GtrModelName::K80, &mut branch_lengths)?;
    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);

    assert_eq!(
      0,
      find_zero_optimal_internal_edges(&graph, &reconstruction, &branch_lengths).len(),
      "precondition: no zero-length internal edges before optimization"
    );

    optimize_edge_lengths(&graph, total_length, &contributions, &indels, method, &mut branch_lengths)?;

    let zero_edges = find_zero_optimal_internal_edges(&graph, &reconstruction, &branch_lengths);
    assert_eq!(
      2,
      zero_edges.len(),
      "{method:?}: expected 2 zero-optimal internal edges after optimization, got {}",
      zero_edges.len()
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rises_to_local_maximum(   0.05, 0.16657011017500076)]
  #[case::descends_to_local_maximum(0.6,  0.16657011017500076)]
  #[case::rises_to_interval_end(    2.0,  3.01               )]
  #[case::reaches_global_maximum(   5.0,  6.369747058293658  )]
  #[trace]
  fn test_dispatch_zero_boundary_solve_newton_never_worse_than_start_and_converges_on_dinh_matsen_k80(
    #[case] t0: f64,
    #[case] expected: f64,
  ) {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let one_mutation = 0.01;

    let solution = solve_newton(t0, &contributions, IndelTerm::new(0, 0.0), 0.0, one_mutation).unwrap();

    assert_eq!(Convergence::Converged, solution.convergence);
    let lh_start = evaluate_partitions_log_lh_only(&contributions, t0).unwrap().value();
    let lh_result = evaluate_partitions_log_lh_only(&contributions, solution.branch_length).unwrap().value();
    assert!(
      lh_result >= lh_start,
      "result t={} has log_lh {lh_result} below the start t0={t0} with {lh_start}",
      solution.branch_length
    );
    pretty_assert_abs_diff_eq!(
      expected,
      solution.branch_length,
      epsilon = 1e-9,
      "expected: stationary points t = -2 ln(x) of P(x) = p1(x) p2(x), x = exp(-t/2), from the roots of the cubic P'(x) of the exact fixture coefficients (numpy): local maximum 0.16657, local minimum 1.2518, global maximum 6.3697; from t0 = 2 the surface rises up to the end 1.5 t0 + one_mutation = 3.01 of the search interval"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_solve_newton_keeps_feasible_start_below_minimum_branch_length() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let (t0, min_branch_length, one_mutation) = (0.1, 0.5, 0.01);

    let solution = solve_newton(
      t0,
      &contributions,
      IndelTerm::new(0, 0.0),
      min_branch_length,
      one_mutation,
    )
    .unwrap();

    assert_eq!(Convergence::Converged, solution.convergence);
    let lh_start = evaluate_partitions_log_lh_only(&contributions, t0).unwrap().value();
    let lh_result = evaluate_partitions_log_lh_only(&contributions, solution.branch_length)
      .unwrap()
      .value();
    assert!(
      lh_result >= lh_start,
      "result t={} has log_lh {lh_result} below the start t0={t0} with {lh_start}",
      solution.branch_length
    );
    pretty_assert_abs_diff_eq!(
      0.16657011017500076,
      solution.branch_length,
      epsilon = 1e-9,
      "expected: local maximum t = 0.16657 of the Dinh and Matsen K80 surface (numpy roots of the cubic P'(x), see the test above), below the minimum branch length 0.5 and above the start 0.1"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_uncertified_returns_grid_point_beating_zero() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let one_mutation = 0.01;
    let start = 0.6;
    let lh_zero = evaluate_partitions_log_lh_only(&contributions, 0.0).unwrap().value();

    let result = reconcile(0.0, start, &contributions, IndelTerm::new(0, 0.0), one_mutation);

    let lh_result = evaluate_partitions_log_lh_only(&contributions, result).unwrap().value();
    assert!(
      lh_result > lh_zero,
      "the grid point must beat zero: log_lh({result})={lh_result} vs log_lh(0)={lh_zero}"
    );
    pretty_assert_ulps_eq!(
      0.16285785765815017,
      result,
      max_ulps = 4,
      "expected the best of 100 log-spaced points on [0.001, 0.91] (index 74), next to the local maximum t = 0.16657 of the exact fixture coefficients (numpy evaluation); the global maximum t = 6.37 lies outside the grid"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_certified_passes_through() {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.0, 1.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr));
    let contributions = [contribution];

    assert!(
      contributions[0].certified_decay_rate().is_some(),
      "precondition: JC69 must have a certified single decay rate"
    );

    let result = reconcile(0.0, 0.1, &contributions, IndelTerm::new(0, 0.0), 0.01);
    assert!(
      result == 0.0,
      "reconcile_zero_boundary must pass exact-zero through for a certified decay rate, got {result}"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_mixed_decay_rates_finds_positive_mode() {
    let contributions = mixed_decay_rate_contributions();
    let result = reconcile(0.0, 0.5, &contributions, IndelTerm::new(0, 0.0), 0.01);
    let lh_zero = evaluate_partitions_log_lh_only(&contributions, 0.0).unwrap().value();
    let lh_result = evaluate_partitions_log_lh_only(&contributions, result).unwrap().value();
    assert!(
      lh_result > lh_zero,
      "contributions with different decay rates must leave zero for a better positive length: log_lh({result})={lh_result} vs log_lh(0)={lh_zero}"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_indels_passes_through() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let result = reconcile(0.0, 0.2, &contributions, IndelTerm::new(1, 44.4), 0.01);
    assert!(
      result == 0.0,
      "reconcile_zero_boundary must pass exact-zero through when indels are present, got {result}"
    );
  }

  fn mixed_decay_rate_contributions() -> [OptimizationContribution; 2] {
    let fast = binary_contribution(5.0, array![[0.96, 1.0]]);
    let slow = binary_contribution(0.5, array![[-0.8, 1.0]]);
    [fast, slow]
  }

  fn binary_contribution(mu: f64, coefficients: Array2<f64>) -> OptimizationContribution {
    let gtr = GTR::builder()
      .n_states(2)
      .mu(Mu::Normalized(mu))
      .pi(array![0.5, 0.5])
      .build()
      .unwrap();
    OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr))
  }

  fn reconcile(
    candidate: f64,
    start: f64,
    contributions: &[OptimizationContribution],
    indels: IndelTerm,
    one_mutation: f64,
  ) -> f64 {
    let zero_length = ZeroLength::of_edge(contributions, indels);
    reconcile_zero_boundary(candidate, start, contributions, indels, zero_length, one_mutation).unwrap()
  }
}
