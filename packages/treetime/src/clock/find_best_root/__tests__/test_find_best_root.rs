#[cfg(test)]
mod tests {
  use crate::clock::find_best_root::params::{
    BranchPointOptimizationParams, BrentParams, GoldenSectionParams, GridSearchParams,
  };
  use crate::o;
  use crate::pretty_assert_abs_diff_eq;
  use crate::test_utils::half_residual_sum_of_squares;
  use eyre::Report;
  use helpers::{
    NEGATIVE_RATE_DATES, OPTIMAL_ROOT_CD_SPLIT, POSITIVE_RATE_DATES, RootSearch, dates_negative_rate,
    dates_positive_rate, divs_on_root_ab_edge, divs_on_root_cd_edge, search_root,
  };
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::assert_error;
  use treetime_utils::least_squares::LineFit;

  #[rustfmt::skip]
  #[rstest]
  #[case::grid(                     BranchPointOptimizationParams::Grid(GridSearchParams::default()),                                                     0.1)]
  #[case::grid_with_params(         BranchPointOptimizationParams::grid_with(GridSearchParams { n_points: 51 }),                                           0.14)]
  #[case::brent(                    BranchPointOptimizationParams::Brent(BrentParams::default()),                                                         OPTIMAL_ROOT_CD_SPLIT)]
  #[case::brent_with_params(        BranchPointOptimizationParams::brent_with(BrentParams { brent_max_iters: 25, brent_tolerance: 1e-8 }),                OPTIMAL_ROOT_CD_SPLIT)]
  #[trace]
  fn test_find_best_root_splits_the_root_to_cd_branch(
    #[case] params: BranchPointOptimizationParams,
    #[case] expected_split: f64,
  ) -> Result<(), Report> {
    let expected_chisq = half_residual_sum_of_squares(&POSITIVE_RATE_DATES, &divs_on_root_cd_edge(expected_split));

    let search = search_root(&dates_positive_rate(), &params, true, None)?;

    assert_eq!(Some((o!("root"), o!("CD"))), search.split_edge);
    pretty_assert_abs_diff_eq!(expected_split, search.split.expect("the root search splits a branch"), epsilon = 1e-7);
    pretty_assert_abs_diff_eq!(expected_chisq, search.result.regression().chisq(), epsilon = 1e-12);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::golden_section(            BranchPointOptimizationParams::GoldenSection(GoldenSectionParams::default()))]
  #[case::golden_section_with_params(BranchPointOptimizationParams::golden_section_with(GoldenSectionParams { golden_max_iters: 25, golden_tolerance: 1e-8 }))]
  #[trace]
  fn test_find_best_root_golden_section_splits_the_root_to_cd_branch(
    #[case] params: BranchPointOptimizationParams,
  ) -> Result<(), Report> {
    let expected_chisq = half_residual_sum_of_squares(&POSITIVE_RATE_DATES, &divs_on_root_cd_edge(OPTIMAL_ROOT_CD_SPLIT));

    let search = search_root(&dates_positive_rate(), &params, true, None)?;

    assert_eq!(Some((o!("root"), o!("CD"))), search.split_edge);
    pretty_assert_abs_diff_eq!(OPTIMAL_ROOT_CD_SPLIT, search.split.expect("the root search splits a branch"), epsilon = 1e-6);
    pretty_assert_abs_diff_eq!(expected_chisq, search.result.regression().chisq(), epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_optimization_methods_improve_on_grid() -> Result<(), Report> {
    let dates = dates_positive_rate();
    let grid = search_root(
      &dates,
      &BranchPointOptimizationParams::Grid(GridSearchParams::default()),
      true,
      None,
    )?;
    let brent = search_root(
      &dates,
      &BranchPointOptimizationParams::Brent(BrentParams::default()),
      true,
      None,
    )?;
    let golden = search_root(
      &dates,
      &BranchPointOptimizationParams::GoldenSection(GoldenSectionParams::default()),
      true,
      None,
    )?;

    let chisq = |search: &RootSearch| search.result.regression().chisq();
    assert!(chisq(&brent) <= chisq(&grid), "Brent should achieve <= chisq than grid");
    assert!(
      chisq(&golden) <= chisq(&grid),
      "Golden section should achieve <= chisq than grid"
    );
    assert_eq!(grid.split_edge, brent.split_edge, "Brent should find same edge as grid");
    assert_eq!(
      grid.split_edge, golden.split_edge,
      "Golden section should find same edge as grid"
    );
    Ok(())
  }

  #[test]
  fn test_find_best_root_fixed_rate_minimizes_the_residuals_around_the_fixed_slope() -> Result<(), Report> {
    let brent = BranchPointOptimizationParams::Brent(BrentParams::default());
    let rate = 1e-3;

    let fixed = search_root(&dates_positive_rate(), &brent, true, Some(rate))?;
    let estimated = search_root(&dates_positive_rate(), &brent, true, None)?;

    assert_eq!(Some((o!("root"), o!("AB"))), fixed.split_edge);
    pretty_assert_abs_diff_eq!(
      0.1675,
      fixed.split.expect("the root search splits a branch"),
      epsilon = 1e-7
    );
    assert_eq!(Some((o!("root"), o!("CD"))), estimated.split_edge);
    let intercept = (0.92 - rate * 8057.0) / 4.0;
    pretty_assert_abs_diff_eq!(
      intercept,
      fixed.result.into_clock_fit()?.model.intercept(),
      epsilon = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_find_best_root_force_positive_true_rejects_negative_rate() -> Result<(), Report> {
    let result = search_root(
      &dates_negative_rate(),
      &BranchPointOptimizationParams::Grid(GridSearchParams::default()),
      true,
      None,
    );

    assert_error!(
      result,
      "Clock rate is negative for all root positions. The data may lack temporal signal. Please specify --clock-rate explicitly."
    );
    Ok(())
  }

  #[test]
  fn test_find_best_root_force_positive_false_accepts_negative_rate() -> Result<(), Report> {
    let expected_split = 0.0125;
    let expected_chisq = half_residual_sum_of_squares(&NEGATIVE_RATE_DATES, &divs_on_root_ab_edge(expected_split));

    let search = search_root(
      &dates_negative_rate(),
      &BranchPointOptimizationParams::Brent(BrentParams::default()),
      false,
      None,
    )?;

    let regression = search.result.regression();
    assert!(
      regression.clock_rate() < 0.0,
      "rate should be negative for this test graph, got {:.6e}",
      regression.clock_rate()
    );
    assert_eq!(Some((o!("root"), o!("AB"))), search.split_edge);
    pretty_assert_abs_diff_eq!(
      expected_split,
      search.split.expect("the root search splits a branch"),
      epsilon = 1e-7
    );
    pretty_assert_abs_diff_eq!(expected_chisq, regression.chisq(), epsilon = 1e-12);
    Ok(())
  }

  #[test]
  fn test_find_best_root_force_positive_accepts_fixed_positive_rate() -> Result<(), Report> {
    let rate = 1e-3;

    let search = search_root(
      &dates_negative_rate(),
      &BranchPointOptimizationParams::Brent(BrentParams::default()),
      true,
      Some(rate),
    )?;

    assert_eq!(Some((o!("root"), o!("AB"))), search.split_edge);
    pretty_assert_abs_diff_eq!(
      0.225,
      search.split.expect("the root search splits a branch"),
      epsilon = 1e-7
    );
    let fit = search.result.into_clock_fit()?;
    let intercept = (0.92 - rate * 8054.0) / 4.0;
    pretty_assert_abs_diff_eq!(intercept, fit.model.intercept(), epsilon = 1e-12);
    let dates: Vec<f64> = fit.points.iter().map(|point| point.date.expect("dated leaf")).collect();
    let divs: Vec<f64> = fit.points.iter().map(|point| point.div).collect();
    let estimated_slope = LineFit::least_squares(&dates, &divs).slope;
    assert!(
      estimated_slope < 0.0,
      "the estimated rate at the fixed-rate root should be negative, got {estimated_slope:.6e}"
    );
    Ok(())
  }

  mod helpers {
    use crate::clock::clock_regression::{
      ClockRerootResult, ClockTree, ClockVarianceParams, estimate_clock_model_with_reroot_policy,
    };
    use crate::clock::clock_state::ClockInputs;
    use crate::clock::find_best_root::params::BranchPointOptimizationParams;
    use crate::clock::reroot::RerootParams;
    use crate::o;
    use crate::progress::NoopProgress;
    use eyre::Report;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_io::nwk::nwk_read;

    pub(super) const LEAF_NAMES: [&str; 4] = ["A", "B", "C", "D"];

    pub(super) const POSITIVE_RATE_DATES: [f64; 4] = [2013.0, 2022.0, 2017.0, 2005.0];

    pub(super) const NEGATIVE_RATE_DATES: [f64; 4] = [2017.0, 2005.0, 2010.0, 2022.0];

    pub(super) const OPTIMAL_ROOT_CD_SPLIT: f64 = 6.18 / 45.0;

    const ROOT_TO_TIP_DIVS: [f64; 4] = [0.2, 0.3, 0.25, 0.17];

    const ROOT_AB_LENGTH: f64 = 0.1;

    const ROOT_CD_LENGTH: f64 = 0.05;

    pub(super) struct RootSearch {
      pub result: ClockRerootResult,
      pub split_edge: Option<(String, String)>,
      pub split: Option<f64>,
    }

    pub(super) fn dates_positive_rate() -> BTreeMap<String, f64> {
      named_dates(POSITIVE_RATE_DATES)
    }

    pub(super) fn dates_negative_rate() -> BTreeMap<String, f64> {
      named_dates(NEGATIVE_RATE_DATES)
    }

    pub(super) fn divs_on_root_cd_edge(split: f64) -> [f64; 4] {
      shifted_divs(ROOT_CD_LENGTH * split)
    }

    pub(super) fn divs_on_root_ab_edge(split: f64) -> [f64; 4] {
      shifted_divs(-ROOT_AB_LENGTH * split)
    }

    pub(super) fn search_root(
      dates: &BTreeMap<String, f64>,
      optimization_params: &BranchPointOptimizationParams,
      force_positive_rate: bool,
      clock_rate: Option<f64>,
    ) -> Result<RootSearch, Report> {
      let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;".as_slice())?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let times = graph
        .get_leaves()
        .map(|node| (node.key(), dates.get(names[&node.key()].as_deref().unwrap()).copied()))
        .collect();
      let edge_names: BTreeMap<_, _> = graph
        .get_edges()
        .map(|edge| {
          let name = |key| names[&key].clone().unwrap_or_else(|| o!("unnamed"));
          (edge.key(), (name(edge.source()), name(edge.target())))
        })
        .collect();
      let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
      let reroot_params = RerootParams {
        force_positive_rate,
        ..RerootParams::default()
      };

      let (_, result) = estimate_clock_model_with_reroot_policy(
        ClockTree {
          graph,
          branch_lengths: nwk_parsed.branch_lengths,
          inputs,
        },
        &BTreeSet::new(),
        &ClockVarianceParams::default(),
        clock_rate,
        false,
        optimization_params,
        &reroot_params,
        None,
        &names,
        &NoopProgress,
      )?;
      let edge_split = result.reroot_result().and_then(|reroot| reroot.edge_split.as_ref());
      let split_edge = edge_split.map(|split| edge_names[&split.old_edge_key].clone());
      let split = edge_split.map(|split| split.split_position);
      Ok(RootSearch {
        result,
        split_edge,
        split,
      })
    }

    fn named_dates(dates: [f64; 4]) -> BTreeMap<String, f64> {
      LEAF_NAMES.iter().map(|name| o!(*name)).zip(dates).collect()
    }

    fn shifted_divs(toward_cd: f64) -> [f64; 4] {
      let [a, b, c, d] = ROOT_TO_TIP_DIVS;
      [a + toward_cd, b + toward_cd, c - toward_cd, d - toward_cd]
    }
  }
}
