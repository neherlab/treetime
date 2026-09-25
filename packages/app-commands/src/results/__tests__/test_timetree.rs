#[cfg(test)]
mod tests {
  use crate::json_float::JsonFloat;
  use crate::results::__tests__::test_tree::tests::helpers::fixture;
  use crate::results::clock::{ClockLine, RootToTip, RootToTipPoint, root_to_tip};
  use crate::results::methods::timetree_methods;
  use crate::results::timetree::{
    CoalescentPrior, RelaxedClock, TimetreeEstimates, TimetreeOutputs, coalescent_prior, timetree_results,
  };
  use crate::results::tree::{DateInterval, ResultTree};
  use eyre::Report;
  use helpers::{clock_model, config, fixed_clock_model, metrics, node_data_clock};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime::clock::rtt::{ClockDateSource, ClockRegressionResult};
  use treetime::o;

  #[test]
  fn test_timetree_estimates_come_from_tree_clock_model_node_data_and_tracelog() -> Result<(), Report> {
    let tree = ResultTree::from_auspice(&fixture())?;
    let model = clock_model();
    let clock = node_data_clock();
    let trace = [metrics(-120.0), metrics(-110.5)];

    let results = timetree_results(
      &TimetreeOutputs {
        tree: Some(&tree),
        clock_model: Some(&model),
        clock_rows: None,
        node_data_clock: Some(&clock),
        trace: &trace,
        coalescent: &[],
      },
      &config(|_| {}),
    );

    let expected = TimetreeEstimates {
      root_date: Some(2010.0),
      root_interval: Some(DateInterval {
        lower: 2009.0,
        upper: 2011.0,
        days: 730.0,
      }),
      root_near_interval_edge: false,
      clock_rate: Some(1e-3),
      clock_rate_std: Some(2e-4),
      clock_rate_fixed: false,
      r: Some(0.9),
      r_squared: Some(0.9 * 0.9),
      samples: 3,
      excluded_samples: 1,
      coalescent_prior: CoalescentPrior::None,
      relaxed_clock: None,
      log_likelihood: Some(JsonFloat(-110.5)),
      iterations: 2,
    };
    assert_eq!(Some(expected), results.estimates);
    assert_eq!(2, results.iterations.len());
    Ok(())
  }

  #[test]
  fn test_timetree_estimates_mark_a_root_date_at_the_interval_edge() -> Result<(), Report> {
    let mut auspice = fixture();
    auspice.tree.node_attrs.num_date.as_mut().unwrap().confidence = Some([2009.95, 2011.95]);
    let tree = ResultTree::from_auspice(&auspice)?;

    let results = timetree_results(
      &TimetreeOutputs {
        tree: Some(&tree),
        clock_model: None,
        clock_rows: None,
        node_data_clock: None,
        trace: &[],
        coalescent: &[],
      },
      &config(|_| {}),
    );

    assert_eq!(
      Some(true),
      results.estimates.map(|estimates| estimates.root_near_interval_edge)
    );
    Ok(())
  }

  #[test]
  fn test_timetree_estimates_take_the_standard_deviation_of_a_fixed_rate_from_the_config() -> Result<(), Report> {
    let tree = ResultTree::from_auspice(&fixture())?;
    let model = fixed_clock_model();
    let clock = node_data_clock();

    let results = timetree_results(
      &TimetreeOutputs {
        tree: Some(&tree),
        clock_model: Some(&model),
        clock_rows: None,
        node_data_clock: Some(&clock),
        trace: &[],
        coalescent: &[],
      },
      &config(|config| config.clock_std_dev = Some(1e-4)),
    );

    let estimates = results.estimates.unwrap();
    assert_eq!(
      (Some(8e-4), Some(1e-4), true, None),
      (
        estimates.clock_rate,
        estimates.clock_rate_std,
        estimates.clock_rate_fixed,
        estimates.r
      )
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::none(     (None,      false, false), CoalescentPrior::None)]
  #[case::fixed(    (Some(0.5), false, false), CoalescentPrior::Fixed { tc: 0.5 })]
  #[case::optimized((Some(0.5), true,  false), CoalescentPrior::Optimized)]
  #[case::skyline(  (None,      true,  true),  CoalescentPrior::Skyline { points: 20, stiffness: 2.0 })]
  #[trace]
  fn test_coalescent_prior_from_config(
    #[case] (tc, optimized, skyline): (Option<f64>, bool, bool),
    #[case] expected: CoalescentPrior,
  ) {
    let config = config(|config| {
      config.coalescent = tc;
      config.coalescent_opt = optimized;
      config.coalescent_skyline = skyline;
      config.skyline_n_points = 20;
      config.skyline_stiffness = 2.0;
    });

    assert_eq!(expected, coalescent_prior(&config));
  }

  #[test]
  fn test_timetree_methods_describe_an_estimated_rate() -> Result<(), Report> {
    let tree = ResultTree::from_auspice(&fixture())?;
    let model = clock_model();
    let clock = node_data_clock();
    let config = config(|_| {});
    let estimates = timetree_results(
      &TimetreeOutputs {
        tree: Some(&tree),
        clock_model: Some(&model),
        clock_rows: None,
        node_data_clock: Some(&clock),
        trace: &[],
        coalescent: &[],
      },
      &config,
    )
    .estimates
    .unwrap();

    let expected = indoc! {"
    A time-scaled phylogeny of 3 samples was inferred with TreeTime 1.2.3 (timetree command). \
    The clock rate was estimated at 1.00e-3 substitutions per site per year (standard deviation 2.00e-4). \
    Samples whose root-to-tip residual exceeded 3 interquartile distances were treated as clock outliers. \
    1 of 3 samples had no usable date or were clock outliers and did not constrain the clock model. \
    The root was dated to 2010-01-01 (90% interval 2009-01-01 to 2011-01-01). \
    Please cite: Sagulenko P, Puller V, Neher RA. TreeTime: Maximum-likelihood phylodynamic analysis. \
    Virus Evolution 4 (2018), vex042."};
    assert_eq!(expected, timetree_methods("1.2.3", &config, &estimates));
    Ok(())
  }

  #[test]
  fn test_timetree_methods_describe_a_fixed_rate_prior_and_relaxed_clock() -> Result<(), Report> {
    let tree = ResultTree::from_auspice(&fixture())?;
    let model = fixed_clock_model();
    let config = config(|config| {
      config.clock_std_dev = Some(1e-4);
      config.clock_filter = 0.0;
      config.coalescent_skyline = true;
      config.skyline_n_points = 20;
      config.skyline_stiffness = 2.0;
      config.relax = vec![1.0, 0.5];
    });
    let mut estimates = timetree_results(
      &TimetreeOutputs {
        tree: Some(&tree),
        clock_model: Some(&model),
        clock_rows: None,
        node_data_clock: None,
        trace: &[],
        coalescent: &[],
      },
      &config,
    )
    .estimates
    .unwrap();
    estimates.root_interval = None;

    let expected = indoc! {"
    A time-scaled phylogeny of 3 samples was inferred with TreeTime 1.2.3 (timetree command). \
    The clock rate was fixed at 8.00e-4 substitutions per site per year with standard deviation 1.00e-4. \
    1 of 3 samples had no usable date or were clock outliers and did not constrain the clock model. \
    A coalescent prior was used (skyline, 20 points, stiffness 2). \
    A relaxed clock was used (slack 1, coupling 0.5). \
    The root was dated to 2010-01-01. \
    Please cite: Sagulenko P, Puller V, Neher RA. TreeTime: Maximum-likelihood phylodynamic analysis. \
    Virus Evolution 4 (2018), vex042."};
    assert_eq!(
      (Some(RelaxedClock {
        slack: 1.0,
        coupling: 0.5
      }),),
      (estimates.relaxed_clock,)
    );
    assert_eq!(expected, timetree_methods("1.2.3", &config, &estimates));
    Ok(())
  }

  #[test]
  fn test_root_to_tip_keeps_the_samples_of_the_tree_with_residuals_in_days() -> Result<(), Report> {
    let tree = ResultTree::from_auspice(&fixture())?;
    let model = clock_model();
    let row = |name: &str, date: Option<f64>, predicted_date: f64, is_outlier: bool| ClockRegressionResult {
      name: Some(name.to_owned()),
      div: 0.01,
      date,
      predicted_date,
      clock_deviation: None,
      is_outlier,
      is_leaf: true,
      date_source: date.map(|_| ClockDateSource::Input),
    };
    let rows = [
      row("AB", None, 2012.0, false),
      row("A", Some(2015.0), 2014.0, false),
      row("C", Some(2014.0), 2015.0, true),
    ];

    let expected = RootToTip {
      points: vec![
        RootToTipPoint {
          name: o!("A"),
          date: Some(2015.0),
          date_source: Some(ClockDateSource::Input),
          div: 0.01,
          predicted_date: 2014.0,
          residual_days: Some(365.0),
          outlier: false,
        },
        RootToTipPoint {
          name: o!("C"),
          date: Some(2014.0),
          date_source: Some(ClockDateSource::Input),
          div: 0.01,
          predicted_date: 2015.0,
          residual_days: Some(-365.0),
          outlier: true,
        },
      ],
      line: Some(ClockLine {
        rate: 1e-3,
        intercept: -2.0,
      }),
    };
    assert_eq!(expected, root_to_tip(Some(&tree), Some(&model), &rows));
    Ok(())
  }

  mod helpers {
    use crate::commands::timetree::args::TreetimeTimetreeArgsRaw;
    use indoc::indoc;
    use std::collections::BTreeMap;
    use treetime::clock::clock_model::ClockModel;
    use treetime::timetree::convergence::metrics::ConvergenceMetrics;
    use treetime_primitives::LogLh;
    use treetime_utils::io::json::json_read_str;
    use util_augur_node_data_json::AugurNodeDataJsonClock;

    pub(super) fn config(change: impl FnOnce(&mut TreetimeTimetreeArgsRaw)) -> TreetimeTimetreeArgsRaw {
      let mut config = TreetimeTimetreeArgsRaw::default();
      change(&mut config);
      config
    }

    pub(super) fn clock_model() -> ClockModel {
      json_read_str(indoc! {r#"{
        "clock_rate": 0.001,
        "intercept": -2.0,
        "stats": {
          "estimated": {
            "chisq": 0.0,
            "r_val": 0.9,
            "hessian": [[1.0, 0.0], [0.0, 1.0]],
            "cov": [[1.0, 0.0], [0.0, 1.0]]
          }
        }
      }"#})
      .unwrap()
    }

    pub(super) fn fixed_clock_model() -> ClockModel {
      json_read_str(r#"{ "clock_rate": 0.0008, "intercept": -1.6, "stats": "fixed" }"#).unwrap()
    }

    pub(super) fn node_data_clock() -> AugurNodeDataJsonClock {
      AugurNodeDataJsonClock {
        rate: 1e-3,
        intercept: -2.0,
        rtt_tmrca: 2000.0,
        cov: None,
        rate_std: Some(2e-4),
        other: BTreeMap::new(),
      }
    }

    pub(super) fn metrics(log_lh_total: f64) -> ConvergenceMetrics {
      ConvergenceMetrics {
        n_diff: 0,
        n_resolved: 0,
        max_time_change: Some(0.1),
        rms_time_change: Some(0.05),
        log_lh_seq: Some(LogLh::new(log_lh_total + 10.0)),
        log_lh_pos: Some(LogLh::new(-10.0)),
        log_lh_coal: None,
        log_lh_total: Some(LogLh::new(log_lh_total)),
      }
    }
  }
}
