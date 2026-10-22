#[cfg(test)]
mod tests {
  use crate::test_utils::PoissonCherry;
  use eyre::Report;
  use helpers::{LAMBDA, TIP_DATE, cherry_input};
  use ndarray::array;
  use rstest::rstest;
  use treetime_utils::cast::ToF64Exact;
  use treetime_utils::{assert_error, pretty_assert_abs_diff_eq};

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_counts_without_prior((0, 0), None     )]
  #[case::zero_counts_with_prior(   (0, 0), Some(1.0))]
  #[case::one_count_without_prior(  (1, 0), None     )]
  #[case::counts_without_prior(     (2, 5), None     )]
  #[case::counts_with_prior(        (3, 1), Some(1.0))]
  #[case::one_zero_count_with_prior((0, 4), Some(0.5))]
  #[trace]
  fn test_reference_log_z_matches_the_poisson_cherry_log_z(
    #[case] counts: (u64, u64),
    #[case] tc: Option<f64>,
  ) -> Result<(), Report> {
    let cherry = PoissonCherry { counts, lambda: LAMBDA, tc };

    let actual = cherry_input(&cherry, &|_| None)?.run()?;

    pretty_assert_abs_diff_eq!(cherry.log_z(), actual.log_z, epsilon = 1e-9);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero_counts_without_prior((0, 0), None     )]
  #[case::zero_counts_with_prior(   (0, 0), Some(1.0))]
  #[case::one_count_without_prior(  (1, 0), None     )]
  #[case::counts_without_prior(     (2, 5), None     )]
  #[case::counts_with_prior(        (3, 1), Some(1.0))]
  #[case::one_zero_count_with_prior((0, 4), Some(0.5))]
  #[trace]
  fn test_reference_log_z_root_age_matches_the_gamma_posterior_of_the_poisson_cherry(
    #[case] counts: (u64, u64),
    #[case] tc: Option<f64>,
  ) -> Result<(), Report> {
    let cherry = PoissonCherry { counts, lambda: LAMBDA, tc };
    let root_age = cherry.root_age();

    let actual = cherry_input(&cherry, &|_| None)?.run()?;

    let mean_age = TIP_DATE - actual.root_time_mean;
    let total_count = (counts.0 + counts.1).to_f64_exact();
    let slope_from_mean_age = total_count - 2.0 * LAMBDA * mean_age;
    pretty_assert_abs_diff_eq!(
      array![root_age.mean(), root_age.variance(), cherry.log_z_slope()],
      array![mean_age, actual.root_time_variance, slope_from_mean_age],
      epsilon = 1e-9
    );
    Ok(())
  }

  #[test]
  fn test_reference_log_z_rejects_a_tip_date_off_the_grid() -> Result<(), Report> {
    let cherry = PoissonCherry {
      counts: (0, 0),
      lambda: LAMBDA,
      tc: None,
    };
    let input = cherry_input(&cherry, &|name| (name == "B").then_some(TIP_DATE - 0.005))?;

    let result = input.run();

    assert_error!(
      result,
      "Tip date 2019.995 lies 0.5000000000109139 grid steps of 0.01 before the youngest tip 2020; the spacing must divide every difference of tip dates"
    );
    Ok(())
  }

  #[test]
  fn test_reference_log_z_rejects_an_interior_jump_that_converges_below_second_order() -> Result<(), Report> {
    let cherry = PoissonCherry {
      counts: (0, 0),
      lambda: LAMBDA,
      tc: None,
    };
    let mut input = cherry_input(&cherry, &|_| None)?;
    input.replace_branch("A", helpers::interior_jump_branch());

    let result = input.run();

    assert_error!(
      result,
      "Reference log Z converges at the observed order NaN below 1.8: -0.11561844638300034, -0.11581140565621695, -0.11571480847133309 at the spacings h, h/2, h/4"
    );
    Ok(())
  }

  mod helpers {
    use crate::test_utils::{
      AnalyticBranch, BranchTime, NodeFactors, PoissonCherry, ReferenceInput, ReferenceLogZ, find_node_key_by_name,
      parent_edge_key, reference_log_z,
    };
    use eyre::Report;
    use std::collections::BTreeMap;
    use std::sync::Arc;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::node_names::NodeNames;
    use treetime_io::nwk::nwk_read;

    pub(super) const LAMBDA: f64 = 1.0;
    pub(super) const TIP_DATE: f64 = 2020.0;
    const SPACING: f64 = 0.01;

    pub(super) struct CherryInput {
      graph: Graph,
      names: NodeNames,
      cherry: PoissonCherry,
      tip_dates: BTreeMap<GraphNodeKey, f64>,
      branches: BTreeMap<GraphEdgeKey, BranchTime>,
    }

    impl CherryInput {
      pub(super) fn run(&self) -> Result<ReferenceLogZ, Report> {
        let root = self.graph.root_key()?;
        let cherry = self.cherry;
        let factor = move |key: GraphNodeKey, _: usize, time: f64| {
          if key == root {
            cherry.root_neg_log_prior(TIP_DATE, time)
          } else {
            0.0
          }
        };
        reference_log_z(&ReferenceInput {
          graph: &self.graph,
          tip_dates: &self.tip_dates,
          branches: &self.branches,
          factors: NodeFactors::Analytic(&factor),
          spacing: SPACING,
        })
      }

      pub(super) fn replace_branch(&mut self, tip: &str, branch: AnalyticBranch) {
        let key = find_node_key_by_name(&self.graph, &self.names, tip).expect("tip must exist");
        self
          .branches
          .insert(parent_edge_key(&self.graph, key), BranchTime::Analytic(branch));
      }
    }

    pub(super) fn cherry_input(
      cherry: &PoissonCherry,
      date_override: &dyn Fn(&str) -> Option<f64>,
    ) -> Result<CherryInput, Report> {
      let parsed = nwk_read(b"(A:1,B:1)root;".as_slice())?;
      let (graph, names) = (parsed.graph, parsed.names);
      let tips = [("A", cherry.counts.0), ("B", cherry.counts.1)];
      let mut tip_dates = BTreeMap::new();
      let mut branches = BTreeMap::new();
      for (name, count) in tips {
        let key = find_node_key_by_name(&graph, &names, name).expect("tip must exist");
        tip_dates.insert(key, date_override(name).unwrap_or(TIP_DATE));
        branches.insert(parent_edge_key(&graph, key), BranchTime::Analytic(cherry.branch(count)));
      }
      Ok(CherryInput {
        graph,
        names,
        cherry: *cherry,
        tip_dates,
        branches,
      })
    }

    pub(super) fn interior_jump_branch() -> AnalyticBranch {
      AnalyticBranch {
        neg_log: Arc::new(|tau: f64| if tau < 0.123_456 { tau } else { tau - 2.0_f64.ln() }),
        monotone_after: 0.0,
      }
    }
  }
}
