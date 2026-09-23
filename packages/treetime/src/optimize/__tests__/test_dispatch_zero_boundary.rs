#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::fitch::create_fitch_partition;
  use crate::ancestral::marginal::branch_lengths_or_zero;
  use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
  use crate::gtr::get_gtr::{GtrModelName, JC69Params, get_gtr_by_name, jc69};
  use crate::optimize::dispatch::{initial_guess_mixed, run_optimize_mixed};
  use crate::optimize::gather::{
    gather_edge_contributions, gather_edge_effective_lengths, gather_edge_indel_counts, gather_edge_sub_counts,
    total_sequence_length,
  };
  use crate::optimize::likelihood::{evaluate_mixed, evaluate_mixed_log_lh_only};
  use crate::optimize::method_newton::{newton_inner, newton_sqrt_inner};
  use crate::optimize::params::BranchOptMethod;
  use crate::optimize::params::ExistingBranchLengths;
  use crate::optimize::run_loop::find_zero_optimal_internal_edges;
  use crate::optimize::run_loop::{marginal_update_dense, marginal_update_sparse};
  use crate::optimize::zero_boundary::{is_zero_branch_optimal, reconcile_zero_boundary};
  use crate::partition::marginal::dense::partition::PartitionMarginalDense;

  use crate::partition::optimize;
  use crate::partition::optimize::contribution::OptimizationContribution;
  use crate::seq::alignment::get_common_length;
  use crate::seq::alignment::node_seq_inputs;
  use eyre::Report;
  use indoc::indoc;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_primitives::AlignmentRecord;

  use ndarray::array;
  use rstest::rstest;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_io::fasta::read_many_fasta_str;
  use treetime_io::nwk::nwk_read_str;

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

  fn setup_identical_partitions(
    graph: &Graph,
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    model: GtrModelName,
    branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
  ) -> Result<(Vec<DenseReconstruction>, Vec<SparseReconstruction>), Report> {
    let aln: Vec<AlignmentRecord> = read_many_fasta_str(IDENTICAL_ALIGNMENT, &Alphabet::default())?
      .into_iter()
      .map(AlignmentRecord::from)
      .collect();

    let dense_partition = PartitionMarginalDense::new(0, Alphabet::new(AlphabetName::Nuc)?, get_common_length(&aln)?);
    let dense_node_states = dense_partition.attach_sequences(graph, &node_seq_inputs(graph, names, aln.clone()))?;
    let dense_partitions = vec![DenseReconstruction::seeded(
      dense_partition,
      get_gtr_by_name(model)?,
      dense_node_states,
    )];

    let fitch = create_fitch_partition(
      graph,
      1,
      Alphabet::new(AlphabetName::Nuc)?,
      &node_seq_inputs(graph, names, aln),
    )?;
    let (sparse_partition, sparse_node_states) = fitch.into_marginal_sparse(graph)?;
    let sparse_partitions = vec![SparseReconstruction::seeded(
      sparse_partition,
      get_gtr_by_name(model)?,
      sparse_node_states,
    )];

    let (dense_partitions, _) =
      marginal_update_dense(graph, &branch_lengths_or_zero(branch_lengths), dense_partitions)?;
    let (sparse_partitions, _) =
      marginal_update_sparse(graph, &branch_lengths_or_zero(branch_lengths), sparse_partitions)?;

    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let indel_counts = gather_edge_indel_counts(graph, &dense_partitions, &sparse_partitions);
    let sub_counts = gather_edge_sub_counts(graph, &dense_partitions, &sparse_partitions)?;
    let effective_lengths = gather_edge_effective_lengths(graph, &dense_partitions, &sparse_partitions)?;
    initial_guess_mixed(
      graph,
      total_length,
      &indel_counts,
      &sub_counts,
      &effective_lengths,
      ExistingBranchLengths::Keep,
      false,
      branch_lengths,
    )?;

    Ok((dense_partitions, sparse_partitions))
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(     BranchOptMethod::Newton)]
  #[case::newton_sqrt(BranchOptMethod::NewtonSqrt)]
  #[case::newton_log( BranchOptMethod::NewtonLog)]
  #[case::brent(      BranchOptMethod::Brent)]
  #[case::brent_sqrt( BranchOptMethod::BrentSqrt)]
  #[case::brent_log(  BranchOptMethod::BrentLog)]
  #[trace]
  fn test_dispatch_zero_boundary_k80_identical_sequences(#[case] method: BranchOptMethod,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(IDENTICAL_TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_identical_partitions(&graph, &names, GtrModelName::K80, &mut branch_lengths)?;
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);

    for (i, edge_ref) in graph.get_edges().enumerate() {
      let bl = branch_lengths[&edge_ref.key()].unwrap();
      assert!(
        bl > 0.0,
        "precondition: edge {i} must start with positive BL, got {bl}"
      );
    }

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;

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
  #[case::tn93( GtrModelName::TN93)]
  #[trace]
  fn test_dispatch_zero_boundary_non_unimodal_models_all_reach_zero(#[case] model: GtrModelName,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(IDENTICAL_TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_identical_partitions(&graph, &names, model, &mut branch_lengths)?;
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);

    assert!(
      !dense_partitions[0].gtr.unimodal_branch_likelihood,
      "precondition: {model:?} must be classified as non-unimodal"
    );

    run_optimize_mixed(
      &graph,
      total_length,
      &contributions,
      &indel_counts,
      BranchOptMethod::BrentSqrt,
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

  #[test]
  fn test_dispatch_zero_boundary_jc69_pre_dispatch_shortcut_reaches_zero() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(IDENTICAL_TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_partitions, sparse_partitions) =
      setup_identical_partitions(&graph, &names, GtrModelName::JC69, &mut branch_lengths)?;
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);

    assert!(
      dense_partitions[0].gtr.unimodal_branch_likelihood,
      "precondition: JC69 must be classified as unimodal"
    );

    let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
    assert!(
      is_zero_branch_optimal(&contributions[&first_edge_key]),
      "precondition: JC69 identical-sequence contributions must trigger the pre-dispatch shortcut"
    );

    run_optimize_mixed(
      &graph,
      total_length,
      &contributions,
      &indel_counts,
      BranchOptMethod::BrentSqrt,
      &mut branch_lengths,
    )?;
    for (i, edge_ref) in graph.get_edges().enumerate() {
      let bl = branch_lengths[&edge_ref.key()].unwrap();
      assert!(
        bl == 0.0,
        "JC69: edge {i} branch length must be exactly 0 after optimization, got {bl}"
      );
    }

    Ok(())
  }

  fn make_dinh_matsen_k80_contribution() -> OptimizationContribution {
    let mut gtr = jc69(JC69Params::default()).unwrap();
    gtr.eigvals = array![-1.0, -0.5, -0.5, 0.0];
    gtr.unimodal_branch_likelihood = false;

    #[rustfmt::skip]
    let coefficients = array![
      [-0.04545042,  0.02261158,  0.02261158, 0.25],
      [ 0.04456328, -0.02228164, -0.02228164, 0.25],
    ];
    OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr))
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_positive_candidate_finds_positive_mode() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let lh_zero = evaluate_mixed_log_lh_only(&contributions, 0.0)
      .expect("valid branch length")
      .value();
    let lh_near_peak = evaluate_mixed_log_lh_only(&contributions, 0.2)
      .expect("valid branch length")
      .value();
    let lh_at_dip = evaluate_mixed_log_lh_only(&contributions, 1.0)
      .expect("valid branch length")
      .value();
    let lh_recovery = evaluate_mixed_log_lh_only(&contributions, 5.0)
      .expect("valid branch length")
      .value();

    assert!(
      lh_near_peak > lh_zero,
      "precondition: log_lh(0.2)={lh_near_peak} > log_lh(0)={lh_zero}"
    );
    assert!(
      lh_near_peak > lh_at_dip,
      "precondition: local max log_lh(0.2)={lh_near_peak} > local dip log_lh(1.0)={lh_at_dip}"
    );
    assert!(
      lh_recovery > lh_at_dip,
      "precondition: recovery log_lh(5.0)={lh_recovery} > local dip log_lh(1.0)={lh_at_dip}"
    );

    let one_mutation = 0.01;
    let candidate = 1.0;
    let lh_candidate = evaluate_mixed_log_lh_only(&contributions, candidate)
      .expect("valid branch length")
      .value();
    assert!(
      lh_candidate < lh_zero,
      "precondition: candidate log_lh({candidate})={lh_candidate} must be worse than log_lh(0)={lh_zero}"
    );

    let branch_length_extent = 0.2;
    let result =
      reconcile_zero_boundary(candidate, branch_length_extent, &contributions, 0, 0.0, one_mutation).unwrap();

    assert!(
      result > 0.0,
      "reconcile_zero_boundary must return a positive mode, not zero, got {result}"
    );
    let lh_result = evaluate_mixed_log_lh_only(&contributions, result)
      .expect("valid branch length")
      .value();
    assert!(
      lh_result > lh_zero,
      "reconciled result must beat zero: log_lh(result)={lh_result} vs log_lh(0)={lh_zero}"
    );
    assert!(
      (0.1..0.4).contains(&result),
      "reconciled result should land near the local max at t ≈ 0.2, got {result}"
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
    let result = reconcile_zero_boundary(candidate, 0.1, &contributions, 0, 0.0, 0.001).unwrap();
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
    let result = reconcile_zero_boundary(candidate, 0.1, &contributions, 0, 0.0, one_mutation).unwrap();
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

    let result =
      reconcile_zero_boundary(candidate, 0.2, &contributions, indel_count, indel_rate, one_mutation).unwrap();
    assert!(
      result.to_bits() == candidate.to_bits(),
      "reconcile_zero_boundary must return candidate unchanged when indel_count > 0, got {result}"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::brent(     BranchOptMethod::Brent)]
  #[case::brent_sqrt(BranchOptMethod::BrentSqrt)]
  #[case::brent_log( BranchOptMethod::BrentLog)]
  #[case::newton_log(BranchOptMethod::NewtonLog)]
  #[trace]
  fn test_dispatch_zero_boundary_topology_cleanup_collects_k80_internal_edges(#[case] method: BranchOptMethod,
  ) -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(IDENTICAL_TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;
    let (dense_partitions, sparse_partitions) = setup_identical_partitions(&graph, &names, GtrModelName::K80, &mut branch_lengths)?;
    let total_length = total_sequence_length(&dense_partitions, &sparse_partitions);
    let contributions = gather_edge_contributions(&graph, &dense_partitions, &sparse_partitions)?;
    let indel_counts = gather_edge_indel_counts(&graph, &dense_partitions, &sparse_partitions);

    assert_eq!(
      0,
      find_zero_optimal_internal_edges(&graph, &sparse_partitions, &branch_lengths).len(),
      "precondition: no zero-length internal edges before optimization"
    );

    run_optimize_mixed(&graph, total_length, &contributions, &indel_counts, method, &mut branch_lengths)?;

    let zero_edges = find_zero_optimal_internal_edges(&graph, &sparse_partitions, &branch_lengths);
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
  #[case::t_0_05(0.05)]
  #[case::t_0_10(0.10)]
  #[case::t_0_15(0.15)]
  #[case::t_0_20(0.20)]
  #[case::t_0_25(0.25)]
  #[case::t_0_30(0.30)]
  #[case::t_0_40(0.40)]
  #[case::t_0_50(0.50)]
  #[case::t_0_60(0.60)]
  #[case::t_0_70(0.70)]
  #[case::t_0_80(0.80)]
  #[case::t_0_90(0.90)]
  #[trace]
  fn test_dispatch_zero_boundary_newton_inner_does_not_clamp_to_zero_on_dinh_matsen_k80(#[case] t0: f64) {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let one_mutation = 0.01;
    let metrics = evaluate_mixed(&contributions, t0).expect("valid branch length");

    let result = newton_inner(t0, &metrics, &contributions, 0, 0.0, 0.0, one_mutation).unwrap();
    assert!(
      result > 0.0,
      "newton_inner from t0={t0} on the Dinh-Matsen K80 surface must return a positive value, got {result}. \
       If this assertion ever fails, `reconcile_zero_boundary` may need an exact-zero entry condition; \
       see the function rustdoc for the full justification."
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_newton_sqrt_inner_clamps_to_zero_on_dinh_matsen_k80() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let one_mutation = 0.01;
    let t0 = 0.6;
    let metrics = evaluate_mixed(&contributions, t0).expect("valid branch length");

    let result = newton_sqrt_inner(t0, &metrics, &contributions, 0, 0.0, 0.0, one_mutation).unwrap();
    assert!(
      result == 0.0,
      "reproduction: newton_sqrt_inner from t0={t0} on the Dinh-Matsen K80 surface must clamp to exactly 0 (else the reconcile zero-candidate gate is no longer necessary), got {result}"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_finds_positive_mode() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let lh_zero = evaluate_mixed_log_lh_only(&contributions, 0.0)
      .expect("valid branch length")
      .value();
    let lh_near_peak = evaluate_mixed_log_lh_only(&contributions, 0.2)
      .expect("valid branch length")
      .value();
    assert!(
      lh_near_peak > lh_zero,
      "precondition: log_lh(0.2)={lh_near_peak} > log_lh(0)={lh_zero}"
    );

    let one_mutation = 0.01;
    let candidate = 0.0;
    let branch_length_extent = 0.6;
    let result =
      reconcile_zero_boundary(candidate, branch_length_extent, &contributions, 0, 0.0, one_mutation).unwrap();

    assert!(
      result > 0.0,
      "reconcile_zero_boundary must reject exact-zero and return a positive mode on a multi-modal surface, got {result}"
    );
    let lh_result = evaluate_mixed_log_lh_only(&contributions, result)
      .expect("valid branch length")
      .value();
    assert!(
      lh_result > lh_zero,
      "reconciled result must beat zero: log_lh(result)={lh_result} vs log_lh(0)={lh_zero}"
    );
    assert!(
      (0.1..0.4).contains(&result),
      "reconciled result should land near the local max at t ≈ 0.2, got {result}"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_unimodal_passes_through() {
    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.0, 1.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize::dense::PartitionContribution::new(coefficients, gtr));
    let contributions = [contribution];

    assert!(
      contributions[0].has_unimodal_branch_likelihood(),
      "precondition: JC69 must be classified as unimodal"
    );

    let result = reconcile_zero_boundary(0.0, 0.1, &contributions, 0, 0.0, 0.01).unwrap();
    assert!(
      result == 0.0,
      "reconcile_zero_boundary must pass exact-zero through for unimodal models, got {result}"
    );
  }

  #[test]
  fn test_dispatch_zero_boundary_reconcile_exact_zero_indels_passes_through() {
    let contributions = [make_dinh_matsen_k80_contribution()];
    let result = reconcile_zero_boundary(0.0, 0.2, &contributions, 1, 44.4, 0.01).unwrap();
    assert!(
      result == 0.0,
      "reconcile_zero_boundary must pass exact-zero through when indels are present, got {result}"
    );
  }
}
