#[cfg(test)]
mod tests {
  use crate::clock::find_best_root::params::{
    BranchPointOptimizationParams, BrentParams, GoldenSectionParams, GridSearchParams,
  };
  use crate::o;
  use crate::pretty_assert_ulps_eq;
  use eyre::Report;
  use helpers::{RootSearch, dates_negative_rate, dates_positive_rate, search_root};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::grid(                     BranchPointOptimizationParams::Grid(GridSearchParams::default()),                                                     0.00026106623586340597)]
  #[case::grid_with_params(         BranchPointOptimizationParams::grid_with(GridSearchParams { n_points: 51 }),                                           0.000_256_025_848_142_593_5)]
  #[case::brent(                    BranchPointOptimizationParams::Brent(BrentParams::default()),                                                         0.000_255_999_999_998_356_5)]
  #[case::brent_with_params(        BranchPointOptimizationParams::brent_with(BrentParams { brent_max_iters: 25, brent_tolerance: 1e-8 }),                0.000_255_999_999_998_356_5)]
  #[case::golden_section(           BranchPointOptimizationParams::GoldenSection(GoldenSectionParams::default()),                                         0.00025599999999690367)]
  #[case::golden_section_with_params(BranchPointOptimizationParams::golden_section_with(GoldenSectionParams { golden_max_iters: 25, golden_tolerance: 1e-8 }), 0.000_255_999_999_998_999_2)]
  #[trace]
  fn test_find_best_root_splits_the_root_to_cd_branch(
    #[case] params: BranchPointOptimizationParams,
    #[case] expected_chisq: f64,
  ) -> Result<(), Report> {
    let search = search_root(&dates_positive_rate(), &params, true, None)?;

    pretty_assert_ulps_eq!(expected_chisq, search.result.regression().chisq(), max_ulps = 4);
    assert_eq!(Some((o!("root"), o!("CD"))), search.split_edge);
    let split = search.split.expect("the root search splits a branch");
    assert!((0.0..=1.0).contains(&split), "split should be in [0, 1]");
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
  fn test_find_best_root_force_positive_true_rejects_negative_rate() -> Result<(), Report> {
    let result = search_root(
      &dates_negative_rate(),
      &BranchPointOptimizationParams::Grid(GridSearchParams::default()),
      true,
      None,
    );

    let err_msg = format!(
      "{:?}",
      result
        .err()
        .expect("force_positive=true should reject all-negative-rate graph")
    );
    assert!(
      err_msg.contains("Clock rate is negative"),
      "Error message should mention negative rate, got: {err_msg}"
    );
    Ok(())
  }

  #[test]
  fn test_find_best_root_force_positive_false_accepts_negative_rate() -> Result<(), Report> {
    let search = search_root(
      &dates_negative_rate(),
      &BranchPointOptimizationParams::Grid(GridSearchParams::default()),
      false,
      None,
    )?;

    let regression = search.result.regression();
    assert!(
      regression.clock_rate() < 0.0,
      "rate should be negative for this test graph, got {:.6e}",
      regression.clock_rate()
    );
    assert!(regression.chisq() >= 0.0, "chisq should be non-negative");
    assert!(regression.chisq().is_finite(), "chisq should be finite");
    if let Some(split) = search.split {
      assert!((0.0..=1.0).contains(&split), "split should be in [0, 1]");
    }
    Ok(())
  }

  #[test]
  fn test_find_best_root_force_positive_accepts_fixed_positive_rate() -> Result<(), Report> {
    let rate = 1e-3;

    let search = search_root(
      &dates_negative_rate(),
      &BranchPointOptimizationParams::Grid(GridSearchParams::default()),
      true,
      Some(rate),
    )?;

    pretty_assert_ulps_eq!(rate, search.result.into_clock_fit()?.model.clock_rate(), max_ulps = 4);
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
    use maplit::btreemap;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_io::nwk::nwk_read_str;

    pub(super) struct RootSearch {
      pub result: ClockRerootResult,
      pub split_edge: Option<(String, String)>,
      pub split: Option<f64>,
    }

    pub(super) fn dates_positive_rate() -> BTreeMap<String, f64> {
      btreemap! {
        o!("A") => 2013.0,
        o!("B") => 2022.0,
        o!("C") => 2017.0,
        o!("D") => 2005.0,
      }
    }

    pub(super) fn dates_negative_rate() -> BTreeMap<String, f64> {
      btreemap! {
        o!("A") => 2017.0,
        o!("B") => 2005.0,
        o!("C") => 2010.0,
        o!("D") => 2022.0,
      }
    }

    pub(super) fn search_root(
      dates: &BTreeMap<String, f64>,
      optimization_params: &BranchPointOptimizationParams,
      force_positive_rate: bool,
      clock_rate: Option<f64>,
    ) -> Result<RootSearch, Report> {
      let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
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
  }
}
