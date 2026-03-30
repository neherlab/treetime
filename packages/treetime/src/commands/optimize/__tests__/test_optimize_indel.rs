#[cfg(test)]
mod tests {
  use crate::commands::optimize::optimize_indel::{estimate_indel_rate, poisson_indel_log_lh};
  use crate::commands::optimize::optimize_unified::{
    initial_guess_mixed, is_zero_branch_optimal, run_optimize_mixed, OptimizationContribution,
  };
  use crate::commands::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    setup_partitions, simple_alignment, TREE_NEWICK,
  };
  use crate::representation::payload::ancestral::GraphAncestral;
  use crate::seq::indel::InDel;
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use treetime_graph::edge::HasBranchLength;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::Seq;

  /// Inject indels onto the first edge in each partition (both dense and sparse).
  fn inject_indels_on_first_edge(
    graph: &GraphAncestral,
    dense_partitions: &[std::sync::Arc<parking_lot::RwLock<crate::representation::partition::marginal_dense::PartitionMarginalDense>>],
    sparse_partitions: &[std::sync::Arc<parking_lot::RwLock<crate::representation::partition::marginal_sparse::PartitionMarginalSparse>>],
    indels: Vec<InDel>,
  ) -> treetime_graph::edge::GraphEdgeKey {
    let first_edge_key = graph.get_edges()[0].read_arc().key();
    for partition in dense_partitions {
      let mut p = partition.write_arc();
      p.edges.get_mut(&first_edge_key).unwrap().indels = indels.clone();
    }
    for partition in sparse_partitions {
      let mut p = partition.write_arc();
      p.edges.get_mut(&first_edge_key).unwrap().indels = indels.clone();
    }
    first_edge_key
  }

  /// `estimate_indel_rate` returns 0 when no edges have indels.
  #[test]
  fn test_optimize_indel_estimate_rate_no_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (_, _, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let rate = estimate_indel_rate(&graph, &mixed_partitions);
    assert_abs_diff_eq!(rate, 0.0, epsilon = 1e-15);
    Ok(())
  }

  /// `estimate_indel_rate` returns total_indels / total_branch_length.
  #[test]
  fn test_optimize_indel_estimate_rate_with_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];
    inject_indels_on_first_edge(&graph, &dense_partitions, &sparse_partitions, indels);

    let rate = estimate_indel_rate(&graph, &mixed_partitions);

    // 2 indels total (1 per partition: dense + sparse), divided by total branch length
    let total_bl: f64 = graph
      .get_edges()
      .iter()
      .map(|e| e.read_arc().payload().read_arc().branch_length().unwrap_or(0.0))
      .sum();
    let expected_rate = 2.0 / total_bl;
    assert_abs_diff_eq!(rate, expected_rate, epsilon = 1e-10);
    Ok(())
  }

  /// When an edge has indels but no substitutions, `initial_guess_mixed` assigns a
  /// non-zero branch length using the Poisson MLE.
  #[test]
  fn test_optimize_indel_initial_guess_nonzero_with_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    let indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?)];
    let first_edge_key = inject_indels_on_first_edge(&graph, &dense_partitions, &sparse_partitions, indels);

    // Force the first edge to have zero substitution-based branch length
    graph.get_edges()[0]
      .write_arc()
      .payload()
      .write_arc()
      .set_branch_length(Some(0.0));

    // Re-estimate indel rate from current branch lengths, then run initial guess
    initial_guess_mixed(&graph, &mixed_partitions)?;

    let bl = graph.get_edges()[0]
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap();

    // If this edge has zero subs but indels, initial_guess should give non-zero
    // (unless the edge also happens to have subs, in which case the sub-based estimate is used)
    let sub_count: usize = mixed_partitions
      .iter()
      .map(|p| p.read_arc().edge_subs(&graph, first_edge_key).map(|s| s.len()))
      .sum::<Result<_, _>>()?;

    if sub_count == 0 {
      assert!(bl > 0.0, "Branch with indels but no subs should have positive initial guess");
    }
    Ok(())
  }

  /// `run_optimize_mixed` assigns non-zero branch length when indels are present
  /// even if substitution evidence alone would favor zero.
  ///
  /// Indels are injected AFTER `update_marginal` because the marginal backward pass
  /// recreates edge partition data from scratch (clearing any pre-existing indels).
  #[test]
  fn test_optimize_indel_run_optimize_nonzero_with_indels() -> Result<(), Report> {
    let graph: GraphAncestral = nwk_read_str(TREE_NEWICK)?;
    let aln = simple_alignment()?;
    let (dense_partitions, sparse_partitions, mixed_partitions) = setup_partitions(&graph, &aln)?;

    // Inject indels AFTER setup (which includes update_marginal) to avoid being wiped
    // by the backward pass that recreates DenseEdgePartition from scratch.
    let indels = vec![InDel::del((0, 3), Seq::try_from_str("ACG")?), InDel::del((5, 8), Seq::try_from_str("ACG")?)];
    inject_indels_on_first_edge(&graph, &dense_partitions, &sparse_partitions, indels);

    run_optimize_mixed(&graph, &mixed_partitions)?;

    // The first edge has indels, so its branch length should be positive
    let bl = graph.get_edges()[0]
      .read_arc()
      .payload()
      .read_arc()
      .branch_length()
      .unwrap();
    assert!(bl > 0.0, "Branch with indels should have positive length after optimization");
    Ok(())
  }

  /// The Poisson indel contribution makes the combined second derivative more negative
  /// (more concave), aiding Newton convergence.
  #[test]
  fn test_optimize_indel_poisson_concavity() {
    let k = 3;
    let mu = 5.0;

    // At several branch lengths, second derivative should always be negative
    for &t in &[0.01, 0.1, 0.5, 1.0, 5.0] {
      let metrics = poisson_indel_log_lh(k, mu, t);
      assert!(
        metrics.second_derivative < 0.0,
        "Poisson log-likelihood should be concave for k>0, but at t={t} got d2={:.6}",
        metrics.second_derivative
      );
    }
  }

  /// The Poisson derivative at the MLE (t = k/mu) is zero.
  #[test]
  fn test_optimize_indel_poisson_mle_derivative_zero() {
    for k in 1..=10 {
      let mu = 3.0;
      let t_mle = k as f64 / mu;
      let metrics = poisson_indel_log_lh(k, mu, t_mle);
      assert_abs_diff_eq!(
        metrics.derivative,
        0.0,
        epsilon = 1e-13
      );
    }
  }

  /// The Poisson log-likelihood at the MLE is the maximum.
  #[test]
  fn test_optimize_indel_poisson_mle_is_maximum() {
    let k = 4;
    let mu = 2.0;
    let t_mle = k as f64 / mu;
    let lh_mle = poisson_indel_log_lh(k, mu, t_mle).log_lh;

    // Check nearby points have lower log-likelihood
    for &delta in &[-0.5, -0.1, -0.01, 0.01, 0.1, 0.5] {
      let t = t_mle + delta;
      if t > 0.0 {
        let lh = poisson_indel_log_lh(k, mu, t).log_lh;
        assert!(lh <= lh_mle + 1e-14, "log-lh at t={t} ({lh}) should be <= log-lh at MLE ({lh_mle})");
      }
    }
  }

  /// Numerical derivative matches analytical derivative.
  #[test]
  fn test_optimize_indel_poisson_numerical_derivative() {
    let k = 3;
    let mu = 5.0;
    let t = 0.3;

    let metrics = poisson_indel_log_lh(k, mu, t);

    // First derivative: central difference with small step
    let h1 = 1e-7;
    let lh_plus = poisson_indel_log_lh(k, mu, t + h1).log_lh;
    let lh_minus = poisson_indel_log_lh(k, mu, t - h1).log_lh;
    let numerical_deriv = (lh_plus - lh_minus) / (2.0 * h1);
    assert_abs_diff_eq!(metrics.derivative, numerical_deriv, epsilon = 1e-6);

    // Second derivative: central difference with larger step to avoid cancellation
    // (rounding error in second difference scales as machine_eps / h^2)
    let h2 = 1e-4;
    let lh_plus2 = poisson_indel_log_lh(k, mu, t + h2).log_lh;
    let lh_center = poisson_indel_log_lh(k, mu, t).log_lh;
    let lh_minus2 = poisson_indel_log_lh(k, mu, t - h2).log_lh;
    let numerical_second = (lh_plus2 - 2.0 * lh_center + lh_minus2) / (h2 * h2);
    assert_abs_diff_eq!(metrics.second_derivative, numerical_second, epsilon = 1e-4);
  }

  /// When `indel_count == 0` and substitution derivative is negative,
  /// `is_zero_branch_optimal` returns true (existing behavior preserved).
  #[test]
  fn test_optimize_indel_zero_branch_no_indels_unchanged() {
    use crate::commands::optimize::optimize_dense;
    use crate::gtr::get_gtr::{JC69Params, jc69};
    use ndarray::array;

    let gtr = jc69(JC69Params::default()).unwrap();
    let coefficients = array![[0.0, 1.0, 0.0, 0.0]];
    let contribution = OptimizationContribution::Dense(optimize_dense::PartitionContribution::new(coefficients, gtr));

    // With zero indels, the existing zero-branch check applies
    assert!(is_zero_branch_optimal(&[contribution]));
  }
}
