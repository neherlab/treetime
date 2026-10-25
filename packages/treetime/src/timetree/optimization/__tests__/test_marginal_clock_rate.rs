#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::DateConstraints;
  use crate::pretty_assert_ulps_eq;
  use crate::test_utils::{find_node_key_by_name, parent_edge_key, unknown_branches};
  use crate::timetree::inference::result::BranchLikelihood;
  use crate::timetree::optimization::marginal_clock_rate::{
    MarginalRateInputs, RateWindow, STEP_FRACTION, branches_at_rate, marginal_log_lh_at_rate, step_marginal_clock_rate,
  };
  use eyre::Report;
  use ndarray::Array1;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::{Distribution, NegLog};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_grid::MaxGridPoints;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;
  use treetime_utils::pretty_assert_abs_diff_eq;

  #[test]
  fn test_marginal_clock_rate_matches_closed_form_curve() -> Result<(), Report> {
    let setup = helpers::cherry()?;
    let inputs = setup.inputs();
    let at = |mu: f64| marginal_log_lh_at_rate(&inputs, mu).map(|log_lh| log_lh.expect("root has a distribution"));
    let curve = |mu: f64| -mu.ln() - (mu - 3.0).powi(2);
    for mu in [2.0, 2.5, 4.0] {
      pretty_assert_abs_diff_eq!(curve(mu) - curve(3.0), at(mu)? - at(3.0)?, epsilon = 1e-13);
    }
    Ok(())
  }

  fn parabola_maximum(reference: f64, h: f64) -> f64 {
    let curve = |x: f64| {
      let mu = reference * x.exp();
      -mu.ln() - (mu - 3.0).powi(2)
    };
    let (below, at, above) = (curve(-h), curve(0.0), curve(h));
    h * (below - above) / (2.0 * (above + below - 2.0 * at))
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::first_step(           (2.7, 1.0,   true ), 0.5,   (helpers::BRANCH_GRID_DX, 1e-14))]
  #[case::later_step(           (2.7, 0.5,   false), 0.25,  (0.0005,                  1e-6 ))]
  #[case::at_an_eighth_of_start((2.8, 0.125, false), 0.125, (helpers::BRANCH_GRID_DX, 1e-6 ))]
  #[trace]
  fn test_marginal_clock_rate_steps_the_step_fraction_toward_a_maximum_between_the_probes_and_halves_the_window(
    #[case] (reference, scale, follow_slope): (f64, f64, bool),
    #[case] next_scale: f64,
    #[case] (branch_grid_dx, tolerance): (f64, f64),
  ) -> Result<(), Report> {
    let setup = helpers::cherry_on(branch_grid_dx)?;
    let window = helpers::window(scale, follow_slope);
    let branches = branches_at_rate(&setup.branches, 1.0, reference)?;
    let step = step_marginal_clock_rate(&setup.inputs_at(&branches, reference, window))?;

    let target = parabola_maximum(reference, window.half_width);
    assert!(target.abs() <= window.half_width);
    pretty_assert_abs_diff_eq!(reference * target.exp(), step.target_rate, epsilon = tolerance);
    pretty_assert_abs_diff_eq!(
      reference * (STEP_FRACTION * target).exp(),
      step.clock_rate,
      epsilon = tolerance
    );
    assert_eq!(helpers::window(next_scale, false), step.window);
    Ok(())
  }

  #[test]
  fn test_marginal_clock_rate_later_step_without_a_maximum_between_the_probes_keeps_the_rate_and_window()
  -> Result<(), Report> {
    let setup = helpers::cherry()?;
    let reference = 2.5;
    let window = helpers::window(0.5, false);
    let branches = branches_at_rate(&setup.branches, 1.0, reference)?;
    let step = step_marginal_clock_rate(&setup.inputs_at(&branches, reference, window))?;
    pretty_assert_ulps_eq!(reference, step.target_rate, max_ulps = 0);
    pretty_assert_ulps_eq!(reference, step.clock_rate, max_ulps = 0);
    assert_eq!(window, step.window);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::convex_to_the_bound(    (1.0, false), 1e-15)]
  #[case::maximum_past_the_probes((2.0, true ), 1e-12)]
  #[trace]
  fn test_marginal_clock_rate_first_step_follows_the_slope_out_of_the_window_and_keeps_its_width(
    #[case] (reference, concave): (f64, bool),
    #[case] tolerance: f64,
  ) -> Result<(), Report> {
    let setup = helpers::cherry()?;
    let window = RateWindow::initial();
    let branches = branches_at_rate(&setup.branches, 1.0, reference)?;
    let step = step_marginal_clock_rate(&setup.inputs_at(&branches, reference, window))?;

    let target = if concave {
      parabola_maximum(reference, window.half_width)
    } else {
      2.0_f64.ln()
    };
    assert!(target > window.half_width && target <= 2.0_f64.ln());
    pretty_assert_abs_diff_eq!(reference * target.exp(), step.target_rate, epsilon = tolerance);
    pretty_assert_abs_diff_eq!(
      reference * (STEP_FRACTION * target).exp(),
      step.clock_rate,
      epsilon = tolerance
    );
    assert_eq!(helpers::window(1.0, false), step.window);
    Ok(())
  }

  #[test]
  fn test_marginal_clock_rate_probe_error_names_the_rate_of_the_probe() -> Result<(), Report> {
    let mut setup = helpers::cherry()?;
    let leaf_a = find_node_key_by_name(&setup.graph, &setup.names, "A").expect("leaf not found");
    setup
      .constraints
      .by_node
      .insert(leaf_a, Some(Arc::new(Distribution::<NegLog>::point(2015.0, f64::NAN))));

    let result = step_marginal_clock_rate(&setup.inputs_at(&setup.branches, 0.5, RateWindow::initial()));

    assert_error!(
      result,
      format!(
        "When computing the marginal likelihood of the node times at clock rate 0.5: When re-peaking the time distribution of node {leaf_a}: Cannot re-peak a distribution point: its amplitude is NaN"
      )
    );
    Ok(())
  }

  mod helpers {
    use super::*;
    use treetime_graph::node_names::NodeNames;
    use treetime_utils::cast::{ToF64Exact, f64_to_usize_saturating};

    const PRECISION: f64 = 4.0;

    pub(super) const BRANCH_GRID_DX: f64 = 0.001;

    pub(super) fn window(scale: f64, follow_slope: bool) -> RateWindow {
      RateWindow {
        half_width: RateWindow::initial().half_width * scale,
        follow_slope,
      }
    }

    pub(super) struct Cherry {
      pub graph: Graph,
      pub constraints: DateConstraints,
      pub bad_branches: BTreeMap<GraphNodeKey, bool>,
      pub branches: BTreeMap<treetime_graph::edge::GraphEdgeKey, BranchLikelihood>,
      pub names: NodeNames,
    }

    impl Cherry {
      pub(super) fn inputs(&self) -> MarginalRateInputs<'_> {
        self.inputs_at(&self.branches, 1.0, RateWindow::initial())
      }

      pub(super) fn inputs_at<'a>(
        &'a self,
        branches: &'a BTreeMap<treetime_graph::edge::GraphEdgeKey, BranchLikelihood>,
        reference_rate: f64,
        window: RateWindow,
      ) -> MarginalRateInputs<'a> {
        MarginalRateInputs {
          graph: &self.graph,
          date_constraints: &self.constraints,
          coalescent: None,
          bad_branches: &self.bad_branches,
          branches,
          reference_rate,
          window,
          max_points: MaxGridPoints::default(),
        }
      }
    }

    pub(super) fn cherry() -> Result<Cherry, Report> {
      cherry_on(BRANCH_GRID_DX)
    }

    pub(super) fn cherry_on(branch_grid_dx: f64) -> Result<Cherry, Report> {
      let parsed = nwk_read(b"(A:8.0,B:5.0)root;".as_slice())?;
      let names = parsed.names.clone();
      let graph = parsed.graph;
      let mut constraints = DateConstraints::default();
      let mut branches = unknown_branches(&graph);
      for (name, date, mean) in [("A", 2015.0, 8.0_f64), ("B", 2014.0, 5.0)] {
        let key = find_node_key_by_name(&graph, &names, name).expect("leaf not found");
        constraints
          .by_node
          .insert(key, Some(Arc::new(Distribution::<NegLog>::point(date, 0.0))));
        let n = f64_to_usize_saturating((2.0 * mean / branch_grid_dx).round()) + 1;
        let t = Array1::from_shape_fn(n, |i| branch_grid_dx * i.to_f64_exact());
        let y = t.mapv(|t| 0.5 * PRECISION * (t - mean).powi(2));
        let branch = BranchLikelihood {
          distribution: Some(Arc::new(Distribution::function(t, y)?)),
          time_length: Some(mean),
        };
        branches.insert(parent_edge_key(&graph, key), branch);
      }
      let bad_branches = graph.get_nodes().map(|node| (node.key(), false)).collect();
      Ok(Cherry {
        graph,
        constraints,
        bad_branches,
        branches,
        names,
      })
    }
  }
}
