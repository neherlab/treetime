#[cfg(test)]
mod tests {
  use crate::clock::clock_model::{ClockModel, ClockRegression};
  use crate::clock::clock_regression::{ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy};
  use crate::clock::clock_set::ClockSet;
  use crate::clock::clock_state::ClockInputs;
  use crate::clock::divergence::root_to_node_divergences;
  use crate::clock::find_best_root::params::BranchPointOptimizationParams;
  use crate::clock::reroot::RerootParams;
  use crate::o;
  use crate::progress::NoopProgress;
  use crate::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};
  use eyre::Report;
  use itertools::Itertools;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::{BTreeMap, BTreeSet};
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  fn compute_naive_rate(dates: &BTreeMap<String, f64>, div: &BTreeMap<String, f64>) -> f64 {
    let t: f64 = dates.values().sum();
    let tsq: f64 = dates.values().map(|&x| x * x).sum();
    let dt: f64 = div.iter().map(|(c, div)| div * dates[c]).sum();
    let d: f64 = div.values().sum();
    (dt * 4.0 - d * t) / (tsq * 4.0 - (t * t))
  }

  #[test]
  fn test_clock_naive_rate() -> Result<(), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    let nwk_parsed = nwk_read_str(TREE_4)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let divergences = root_to_node_divergences(&graph, |edge_key| branch_lengths[&edge_key].unwrap_or_default())?;
    let divs: BTreeMap<String, f64> = graph
      .get_leaves()
      .map(|leaf| (names[&leaf.key()].clone().unwrap(), divergences[&leaf.key()]))
      .collect();
    let naive_rate = compute_naive_rate(&dates, &divs);

    let regression = helpers::root_regression(TREE_4, &dates, &ClockVarianceParams::default())?;
    pretty_assert_abs_diff_eq!(naive_rate, regression.clock_rate(), epsilon = 1e-10);

    let options = &ClockVarianceParams {
      variance_factor: 1.0,
      variance_offset: 0.0,
      variance_offset_leaf: 1.0,
    };

    let regression = helpers::root_regression(TREE_4, &dates, options)?;
    pretty_assert_ulps_eq!(0.007710610618916924, regression.clock_rate(), max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_clock_regression_dateless_leaf_does_not_change_root_statistics() -> Result<(), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    let options = ClockVarianceParams::default();
    let expected = helpers::root_regression("(A:0.1,B:0.2,C:0.2,D:0.12)root;", &dates, &options)?;
    let actual = helpers::root_regression("(A:0.1,B:0.2,C:0.2,D:0.12,E:10.0)root;", &dates, &options)?;

    assert_eq!(helpers::regression_bits(&expected), helpers::regression_bits(&actual));
    Ok(())
  }

  #[test]
  fn test_clock_regression_points_are_root_to_tip_distances_of_the_leaves() -> Result<(), Report> {
    let (_, points) = helpers::fit(TREE_4, &dates_4(), &[], true, None)?;

    #[rustfmt::skip]
    let expected = [
      (o!("A"), Some(2013.0), 0.2,  false),
      (o!("B"), Some(2022.0), 0.3,  false),
      (o!("C"), Some(2017.0), 0.25, false),
      (o!("D"), Some(2005.0), 0.17, false),
    ];
    assert_eq!(expected.len(), points.len());
    for ((name, date, div, is_outlier), actual) in expected.iter().zip_eq(&points) {
      assert_eq!((name, date, is_outlier), (&actual.0, &actual.1, &actual.3));
      pretty_assert_ulps_eq!(*div, actual.2, max_ulps = 4);
    }
    Ok(())
  }

  #[test]
  fn test_clock_regression_points_use_clock_lengths_when_the_previous_rate_is_given() -> Result<(), Report> {
    let (_, points) = helpers::fit(TREE_4, &dates_4(), &[], true, Some(0.01))?;

    let (_, input_lengths) = helpers::fit(TREE_4, &dates_4(), &[], true, None)?;

    for ((_, _, div, _), (_, _, input_div, _)) in points.iter().zip_eq(&input_lengths) {
      pretty_assert_ulps_eq!(3.0 * input_div, *div, max_ulps = 4);
    }
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::input_root(               (true,  None))]
  #[case::input_root_clock_lengths( (true,  Some(0.004)))]
  #[case::rerooted(                 (false, None))]
  #[trace]
  fn test_clock_regression_points_refit_to_the_fitted_model(
    #[case] (keep_root, prev_clock_rate): (bool, Option<f64>),
  ) -> Result<(), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
      o!("E") => 2020.0,
    };
    let (model, points) = helpers::fit(TREE_5, &dates, &[], keep_root, prev_clock_rate)?;

    let mut refit = ClockSet::default();
    for (_, date, div, _) in &points {
      refit += ClockSet::leaf_contribution_to_parent(*date, *div, 1.0);
    }
    let refit = ClockRegression::try_from(&refit)?;

    pretty_assert_abs_diff_eq!(model.clock_rate(), refit.clock_rate(), epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(model.intercept(), refit.intercept(), epsilon = 1e-9);
    Ok(())
  }

  #[test]
  fn test_clock_regression_points_flag_the_outliers() -> Result<(), Report> {
    let (_, points) = helpers::fit(TREE_5, &dates_5(), &["E"], true, None)?;

    let flagged = points
      .iter()
      .filter(|(_, _, _, is_outlier)| *is_outlier)
      .map(|(name, ..)| name.as_str())
      .collect_vec();
    assert_eq!(vec!["E"], flagged);
    Ok(())
  }

  #[test]
  #[ignore = "outlier leaves still enter the regression through their parent edge: kb/issues/H-clock-regression-includes-filtered-outliers.md"]
  fn test_clock_regression_outliers_are_left_out_of_the_fit() -> Result<(), Report> {
    let (model, _) = helpers::fit(TREE_5, &dates_5(), &["E"], true, None)?;
    let (model_without_e_date, _) = helpers::fit(TREE_5, &dates_4(), &[], true, None)?;

    pretty_assert_abs_diff_eq!(model_without_e_date.clock_rate(), model.clock_rate(), epsilon = 1e-12);
    Ok(())
  }

  const TREE_4: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  const TREE_5: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,(D:0.12,E:0.3)DE:0.07)CDE:0.05)root:0.01;";

  fn dates_4() -> BTreeMap<String, f64> {
    btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    }
  }

  fn dates_5() -> BTreeMap<String, f64> {
    btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
      o!("E") => 2020.0,
    }
  }

  mod helpers {
    use super::*;

    pub(super) fn fit(
      tree: &str,
      dates: &BTreeMap<String, f64>,
      outlier_names: &[&str],
      keep_root: bool,
      prev_clock_rate: Option<f64>,
    ) -> Result<(ClockModel, Vec<(String, Option<f64>, f64, bool)>), Report> {
      let nwk_parsed = nwk_read_str(tree)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let times = leaf_times(&names, &graph, dates);
      let edge_inputs = branch_lengths
        .iter()
        .map(|(&key, length)| (key, (length.map(|length| length * 100.0), 3.0)))
        .collect();
      let inputs = ClockInputs::from_times(&graph, &times, &edge_inputs);
      let outliers: BTreeSet<GraphNodeKey> = names
        .iter()
        .filter(|(_, name)| name.as_deref().is_some_and(|name| outlier_names.contains(&name)))
        .map(|(key, _)| *key)
        .collect();
      let (_, result) = estimate_clock_model_with_reroot_policy(
        ClockTree {
          graph,
          branch_lengths,
          inputs,
        },
        &outliers,
        &ClockVarianceParams::default(),
        None,
        keep_root,
        &BranchPointOptimizationParams::default(),
        &RerootParams::default(),
        prev_clock_rate,
        &names,
        &NoopProgress,
      )?;
      let fit = result.into_clock_fit()?;
      let points = fit
        .points
        .iter()
        .map(|point| {
          let name = names[&point.key].clone().unwrap();
          (name, point.date, point.div, point.is_outlier)
        })
        .sorted_by(|a, b| a.0.cmp(&b.0))
        .collect_vec();
      Ok((fit.model, points))
    }

    pub(super) fn leaf_times(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      dates: &BTreeMap<String, f64>,
    ) -> BTreeMap<GraphNodeKey, Option<f64>> {
      graph
        .get_leaves()
        .map(|node| {
          let name = names[&node.key()].clone().unwrap();
          (node.key(), dates.get(&name).copied())
        })
        .collect()
    }

    pub(super) fn root_regression(
      tree: &str,
      dates: &BTreeMap<String, f64>,
      options: &ClockVarianceParams,
    ) -> Result<ClockRegression, Report> {
      let nwk_parsed = nwk_read_str(tree)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let times = leaf_times(&names, &graph, dates);
      let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
      let (_, result) = estimate_clock_model_with_reroot_policy(
        ClockTree {
          graph,
          branch_lengths: nwk_parsed.branch_lengths,
          inputs,
        },
        &BTreeSet::new(),
        options,
        None,
        true,
        &BranchPointOptimizationParams::default(),
        &RerootParams::default(),
        None,
        &names,
        &NoopProgress,
      )?;
      Ok(result.regression().clone())
    }

    pub(super) fn regression_bits(regression: &ClockRegression) -> [u64; 4] {
      [
        regression.clock_rate().to_bits(),
        regression.intercept().to_bits(),
        regression.chisq().to_bits(),
        regression.r_squared().to_bits(),
      ]
    }
  }
}
