#[cfg(test)]
pub(super) mod tests {
  use crate::optimize::__tests__::test_convergence::test_convergence_support::tests::{
    TREE_NEWICK, setup_reconstruction, simple_alignment,
  };
  use crate::optimize::dispatch::optimize_edge_lengths;
  use crate::optimize::gather::gather_edge_contributions;
  use crate::optimize::indel::IndelPolicy;
  use crate::optimize::likelihood::{OptimizationMetrics, evaluate_partitions, evaluate_partitions_log_lh_only};
  use crate::optimize::method_newton::{chain_rule_sqrt, newton_tolerance_sqrt};
  use crate::optimize::params::BranchOptMethod;
  use crate::partition::marginal::reconstruction::MarginalReconstruction;
  use crate::seq::indel::InDel;
  use crate::test_utils::{deletion, sparse_reconstruction_mut};
  use eyre::Report;
  use helpers::*;
  use proptest::prelude::*;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read;
  use treetime_primitives::Seq;

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(BranchOptMethod::Newton)]
  #[case::brent( BranchOptMethod::Brent)]
  #[trace]
  fn test_optimize_method_equivalence_no_indels(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let mut branch_lengths = nwk_parsed.branch_lengths;
    let aln = simple_alignment()?;
    let reconstruction = setup_reconstruction(&graph, &names, &aln, &mut branch_lengths)?;
    let total_length = reconstruction.sequence_length();
    let contributions = gather_edge_contributions(&graph, &reconstruction)?;
    let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);

    optimize_edge_lengths(&graph, total_length, &contributions, &indels, method, &mut branch_lengths)?;

    for (i, edge_ref) in graph.get_edges().enumerate() {
      let bl = branch_lengths[&edge_ref.key()].unwrap_or(f64::NAN);
      assert!(bl.is_finite(), "Edge {i}: branch length is not finite ({bl})");
      assert!(bl >= 0.0, "Edge {i}: branch length is negative ({bl})");
    }

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::newton(BranchOptMethod::Newton)]
  #[case::brent( BranchOptMethod::Brent)]
  #[trace]
  fn test_optimize_method_local_optimality(#[case] method: BranchOptMethod) -> Result<(), Report> {
    let OptimizedFirstEdge { graph, reconstruction, indel_rate, branch_length: bl } = optimize_first_edge_with_indels(method, 4)?;
    assert!(bl > 0.0 && bl.is_finite(), "Optimized BL must be positive and finite, got {bl}");

    let lh_opt = eval_combined_first_edge(&graph, &reconstruction, indel_rate, bl)?;

    for &frac in &[0.001, 0.01, 0.1] {
      let delta = bl * frac;
      if bl - delta > 0.0 {
        let lh_below = eval_combined_first_edge(&graph, &reconstruction, indel_rate, bl - delta)?;
        assert!(
          lh_opt >= lh_below - 1e-10,
          "{method:?}: lh at t*={bl} ({lh_opt}) < lh at t*-{frac}*t ({lh_below})"
        );
      }
      let lh_above = eval_combined_first_edge(&graph, &reconstruction, indel_rate, bl + delta)?;
      assert!(
        lh_opt >= lh_above - 1e-10,
        "{method:?}: lh at t*={bl} ({lh_opt}) < lh at t*+{frac}*t ({lh_above})"
      );
    }

    Ok(())
  }

  #[test]
  fn test_optimize_method_newton_sqrt_stationarity() -> Result<(), Report> {
    let OptimizedFirstEdge {
      graph,
      reconstruction,
      indel_rate,
      branch_length: bl,
    } = optimize_first_edge_with_indels(BranchOptMethod::Newton, 2)?;
    let metrics = eval_metrics_first_edge(&graph, &reconstruction, indel_rate, bl)?;
    let s = bl.sqrt();
    let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, metrics.derivative, metrics.second_derivative);

    assert!(
      d2l_ds2 < 0.0,
      "objective must be concave in s at the optimum, d2l_ds2={d2l_ds2}, t*={bl}"
    );
    let implied_step = (dl_ds / d2l_ds2).abs();
    let tol = newton_tolerance_sqrt(s);
    assert!(
      implied_step < tol * 10.0,
      "implied Newton step in s ({implied_step}) exceeds 10x tolerance ({tol}), \
       dl_ds={dl_ds}, d2l_ds2={d2l_ds2}, t*={bl}"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::k1(1)]
  #[case::k2(2)]
  #[case::k4(4)]
  #[trace]
  fn test_optimize_method_cross_method_lh_agreement(#[case] n_indels: usize) -> Result<(), Report> {
    let lh_brent = optimize_first_edge_with_indels(BranchOptMethod::Brent, n_indels)?.log_lh()?;
    let lh_newton = optimize_first_edge_with_indels(BranchOptMethod::Newton, n_indels)?.log_lh()?;

    let lh_diff = (lh_brent - lh_newton).abs();
    assert!(
      lh_diff < 1e-9,
      "Brent lh ({lh_brent}) and Newton lh ({lh_newton}) differ by {lh_diff}"
    );

    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::brent_k1( BranchOptMethod::Brent,  1)]
  #[case::brent_k2( BranchOptMethod::Brent,  2)]
  #[case::brent_k4( BranchOptMethod::Brent,  4)]
  #[case::newton_k1(BranchOptMethod::Newton, 1)]
  #[case::newton_k2(BranchOptMethod::Newton, 2)]
  #[case::newton_k4(BranchOptMethod::Newton, 4)]
  #[trace]
  fn test_optimize_method_positive_with_indels(#[case] method: BranchOptMethod, #[case] n_indels: usize) -> Result<(), Report> {
    let bl = optimize_first_edge_with_indels(method, n_indels)?.branch_length;
    assert!(bl > 0.0, "{method:?} BL with {n_indels} indels must be positive, got {bl}");
    assert!(bl.is_finite(), "{method:?} BL must be finite, got {bl}");
    Ok(())
  }

  mod generators {
    use proptest::prelude::*;
    pub(crate) fn gen_s() -> impl Strategy<Value = f64> {
      1e-6_f64..1e3_f64
    }
    pub(crate) fn gen_dl_dt() -> impl Strategy<Value = f64> {
      -1e6_f64..1e6_f64
    }
    pub(crate) fn gen_d2l_dt2() -> impl Strategy<Value = f64> {
      -1e8_f64..1e8_f64
    }
    pub(crate) fn gen_scalar() -> impl Strategy<Value = f64> {
      prop_oneof![-1e3_f64..-1e-3_f64, 1e-3_f64..1e3_f64]
    }
  }

  pub(crate) mod helpers {
    use super::*;
    use crate::optimize::gather::gather_edge_indel_counts;
    use crate::optimize::indel::poisson_indel_log_lh;

    pub(crate) fn setup_with_indels(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      branch_lengths: &mut BTreeMap<GraphEdgeKey, Option<f64>>,
      n_indels: usize,
    ) -> Result<(MarginalReconstruction, f64), Report> {
      let aln = simple_alignment()?;
      let mut reconstruction = setup_reconstruction(graph, names, &aln, branch_lengths)?;

      let first_edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
      let indels: Vec<InDel> = (0..n_indels)
        .map(|i| deletion((i * 3, i * 3 + 3), Seq::try_from_str("ACG").unwrap()))
        .collect();

      sparse_reconstruction_mut(&mut reconstruction)
        .partition
        .obs_edges
        .get_mut(&first_edge_key)
        .unwrap()
        .indels = indels;

      let indel_rate = IndelPolicy::Include
        .edge_indels(graph, &reconstruction, branch_lengths)
        .rate();

      Ok((reconstruction, indel_rate))
    }

    pub(crate) fn optimize_first_edge_with_indels(
      method: BranchOptMethod,
      n_indels: usize,
    ) -> Result<OptimizedFirstEdge, Report> {
      let nwk_parsed = nwk_read(TREE_NEWICK.as_bytes())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let mut branch_lengths = nwk_parsed.branch_lengths;
      let (reconstruction, indel_rate) = setup_with_indels(&graph, &names, &mut branch_lengths, n_indels)?;
      let total_length = reconstruction.sequence_length();
      let contributions = gather_edge_contributions(&graph, &reconstruction)?;
      let indels = IndelPolicy::Include.edge_indels(&graph, &reconstruction, &branch_lengths);
      optimize_edge_lengths(
        &graph,
        total_length,
        &contributions,
        &indels,
        method,
        &mut branch_lengths,
      )?;
      let branch_length = first_edge_bl(&graph, &branch_lengths);
      Ok(OptimizedFirstEdge {
        graph,
        reconstruction,
        indel_rate,
        branch_length,
      })
    }

    pub(crate) struct OptimizedFirstEdge {
      pub graph: Graph,
      pub reconstruction: MarginalReconstruction,
      pub indel_rate: f64,
      pub branch_length: f64,
    }

    impl OptimizedFirstEdge {
      pub(crate) fn log_lh(&self) -> Result<f64, Report> {
        eval_combined_first_edge(&self.graph, &self.reconstruction, self.indel_rate, self.branch_length)
      }
    }

    pub(crate) fn eval_combined_first_edge(
      graph: &Graph,
      reconstruction: &MarginalReconstruction,
      indel_rate: f64,
      t: f64,
    ) -> Result<f64, Report> {
      let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
      let contributions = gather_edge_contributions(graph, reconstruction)?;
      let indel_counts = gather_edge_indel_counts(graph, reconstruction);
      let edge_contributions = &contributions[&edge_key];
      let indel_count: usize = indel_counts[&edge_key];

      let sub_lh = evaluate_partitions_log_lh_only(edge_contributions, t)
        .expect("valid branch length")
        .value();
      let indel_lh = poisson_indel_log_lh(indel_count, indel_rate, t)
        .expect("valid Poisson parameters")
        .log_lh
        .value();
      Ok(sub_lh + indel_lh)
    }

    pub(crate) fn eval_metrics_first_edge(
      graph: &Graph,
      reconstruction: &MarginalReconstruction,
      indel_rate: f64,
      t: f64,
    ) -> Result<OptimizationMetrics, Report> {
      let edge_key = graph.get_edges().collect::<Vec<_>>()[0].key();
      let contributions = gather_edge_contributions(graph, reconstruction)?;
      let indel_counts = gather_edge_indel_counts(graph, reconstruction);
      let edge_contributions = &contributions[&edge_key];
      let indel_count: usize = indel_counts[&edge_key];

      let mut metrics = evaluate_partitions(edge_contributions, t).expect("valid branch length");
      metrics.add(&poisson_indel_log_lh(indel_count, indel_rate, t).expect("valid Poisson parameters"));
      Ok(metrics)
    }

    pub(crate) fn first_edge_bl(graph: &Graph, branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>) -> f64 {
      branch_lengths[&graph.get_edges().collect::<Vec<_>>()[0].key()].unwrap()
    }

    proptest! {
      #[test]
      fn test_prop_optimize_method_chain_rule_sqrt_formula(
        s in generators::gen_s(),
        dl_dt in generators::gen_dl_dt(),
        d2l_dt2 in generators::gen_d2l_dt2(),
      ) {
        let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, dl_dt, d2l_dt2);
        let expected_dl_ds = 2.0 * s * dl_dt;
        let expected_d2l_ds2 = 4.0 * s * s * d2l_dt2 + 2.0 * dl_dt;
        let dl_tol = 1e-9 * expected_dl_ds.abs().max(1.0);
        let d2l_tol = 1e-9 * expected_d2l_ds2.abs().max(1.0);
        prop_assert!(
          (dl_ds - expected_dl_ds).abs() <= dl_tol,
          "dl_ds: got {dl_ds}, expected {expected_dl_ds}"
        );
        prop_assert!(
          (d2l_ds2 - expected_d2l_ds2).abs() <= d2l_tol,
          "d2l_ds2: got {d2l_ds2}, expected {expected_d2l_ds2}"
        );
      }

      #[test]
      fn test_prop_optimize_method_chain_rule_sqrt_linear(
        s in generators::gen_s(),
        dl_dt in generators::gen_dl_dt(),
        d2l_dt2 in generators::gen_d2l_dt2(),
        k in generators::gen_scalar(),
      ) {
        let (dl_ds, d2l_ds2) = chain_rule_sqrt(s, dl_dt, d2l_dt2);
        let (dl_ds_k, d2l_ds2_k) = chain_rule_sqrt(s, k * dl_dt, k * d2l_dt2);
        let dl_tol = 1e-9 * (k * dl_ds).abs().max(1.0);
        let d2l_tol = 1e-9 * (k * d2l_ds2).abs().max(1.0);
        prop_assert!(
          (dl_ds_k - k * dl_ds).abs() <= dl_tol,
          "linearity violated for dl_ds: {dl_ds_k} vs k*{dl_ds}"
        );
        prop_assert!(
          (d2l_ds2_k - k * d2l_ds2).abs() <= d2l_tol,
          "linearity violated for d2l_ds2: {d2l_ds2_k} vs k*{d2l_ds2}"
        );
      }
    }
  }
}
