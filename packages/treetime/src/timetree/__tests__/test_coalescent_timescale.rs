#[cfg(test)]
mod tests {
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
  use crate::timetree::coalescent::{
    CoalescentBand, CoalescentInputs, CoalescentOutput, CoalescentOutputMode, CoalescentSolve,
  };
  use crate::timetree::coalescent_timescale::{
    CoalescentMode, CoalescentReportBand, CoalescentTcReport, CoalescentTimescale, build_coalescent_output,
    coalescent_mode, estimate_coalescent_tc,
  };
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use maplit::btreemap;
  use ndarray::array;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;
  use treetime_grid::piecewise_constant_fn::PiecewiseConstantFn;
  use treetime_io::dates_csv::{DateConstraint, DatesMap};
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::{o, pretty_assert_array_eq};

  const GEN_PER_YEAR: f64 = 50.0;
  const N_STD: f64 = 2.0;

  #[rustfmt::skip]
  #[rstest]
  #[case::disabled(       None,       false, false, CoalescentMode::Disabled)]
  #[case::fixed(          Some(0.25), false, false, CoalescentMode::Fixed(0.25))]
  #[case::opt_default(    None,       true,  false, CoalescentMode::Constant)]
  #[case::opt_value_ignored(Some(0.25), true, false, CoalescentMode::Constant)]
  #[case::skyline_default(None,       false, true,  CoalescentMode::Skyline)]
  #[case::skyline_over_opt(Some(0.25), true, true,  CoalescentMode::Skyline)]
  #[trace]
  fn test_pipeline_coalescent_mode(
    #[case] coalescent: Option<f64>,
    #[case] coalescent_opt: bool,
    #[case] coalescent_skyline: bool,
    #[case] expected: CoalescentMode,
  ) {
    let actual = coalescent_mode(coalescent, coalescent_opt, coalescent_skyline);

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_pipeline_build_coalescent_output_disabled_returns_none() -> Result<(), Report> {
    let timescale = CoalescentTimescale::constant(1.0);
    let params = SkylineParams {
      n_std: N_STD,
      ..SkylineParams::default()
    };

    let actual = build_coalescent_output(CoalescentMode::Disabled, &timescale, GEN_PER_YEAR, &params)?;

    assert_eq!(None, actual);
    Ok(())
  }

  #[test]
  fn test_pipeline_build_coalescent_output_fixed_emits_one_segment_no_band() -> Result<(), Report> {
    let (graph, constraints) = dated_tree()?;
    let params = SkylineParams {
      n_std: N_STD,
      ..SkylineParams::default()
    };
    let node_times = TimetreeState::seed_from_values(&graph, &constraints).coalescent_node_times()?;
    let timescale = estimate_coalescent_tc(CoalescentMode::Fixed(2.5), &graph, &params, &node_times)?
      .expect("a fixed Tc yields a coalescent timescale");

    let actual = build_coalescent_output(CoalescentMode::Fixed(2.5), &timescale, GEN_PER_YEAR, &params)?
      .expect("a fixed Tc writes a coalescent output");

    let expected = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Fixed,
        n_points: None,
        stiffness: None,
        confidence_n_std: None,
        gen_per_year: GEN_PER_YEAR,
      },
      &CoalescentSolve {
        segment_boundaries: &[2000.0, 2010.0],
        tc_values: &[2.5],
        band: None,
        log_likelihood: None,
      },
    )?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_pipeline_build_coalescent_output_constant_carries_band() -> Result<(), Report> {
    let params = SkylineParams {
      n_std: N_STD,
      ..SkylineParams::default()
    };
    let timescale = CoalescentTimescale {
      distribution: Distribution::constant(3.0),
      schedule: PiecewiseConstantFn::new(array![], array![3.0]),
      report: Some(CoalescentTcReport {
        segment_boundaries: array![2000.0, 2020.0],
        band: Some(CoalescentReportBand {
          lower: array![2.0],
          upper: array![4.0],
        }),
        log_likelihood: Some(-7.5),
      }),
    };
    let actual = build_coalescent_output(CoalescentMode::Constant, &timescale, GEN_PER_YEAR, &params)?;

    let expected = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Constant,
        n_points: None,
        stiffness: None,
        confidence_n_std: Some(N_STD),
        gen_per_year: GEN_PER_YEAR,
      },
      &CoalescentSolve {
        segment_boundaries: &[2000.0, 2020.0],
        tc_values: &[3.0],
        band: Some(CoalescentBand {
          lower: &[2.0],
          upper: &[4.0],
        }),
        log_likelihood: Some(-7.5),
      },
    )?;
    assert_eq!(Some(expected), actual);
    Ok(())
  }

  #[test]
  fn test_pipeline_build_coalescent_output_skyline_multi_segment_band() -> Result<(), Report> {
    let params = SkylineParams {
      n_points: 2,
      stiffness: 3.0,
      n_std: N_STD,
      ..SkylineParams::default()
    };
    let timescale = CoalescentTimescale {
      distribution: Distribution::constant(3.0),
      schedule: PiecewiseConstantFn::new(array![2010.0], array![3.0, 5.0]),
      report: Some(CoalescentTcReport {
        segment_boundaries: array![2000.0, 2010.0, 2020.0],
        band: Some(CoalescentReportBand {
          lower: array![2.0, 4.0],
          upper: array![4.0, 6.0],
        }),
        log_likelihood: Some(-9.0),
      }),
    };
    let actual = build_coalescent_output(CoalescentMode::Skyline, &timescale, GEN_PER_YEAR, &params)?;

    let expected = CoalescentOutput::new(
      CoalescentInputs {
        mode: CoalescentOutputMode::Skyline,
        n_points: Some(2),
        stiffness: Some(3.0),
        confidence_n_std: Some(N_STD),
        gen_per_year: GEN_PER_YEAR,
      },
      &CoalescentSolve {
        segment_boundaries: &[2000.0, 2010.0, 2020.0],
        tc_values: &[3.0, 5.0],
        band: Some(CoalescentBand {
          lower: &[2.0, 4.0],
          upper: &[4.0, 6.0],
        }),
        log_likelihood: Some(-9.0),
      },
    )?;
    assert_eq!(Some(expected), actual);
    Ok(())
  }

  fn dated_tree() -> Result<(Graph, DateConstraints), Report> {
    let dates: DatesMap = btreemap! {
      o!("root") => Some(DateConstraint::exact(2000.0)),
      o!("x")    => Some(DateConstraint::exact(2005.0)),
      o!("a")    => Some(DateConstraint::exact(2010.0)),
      o!("b")    => Some(DateConstraint::exact(2010.0)),
      o!("c")    => Some(DateConstraint::exact(2010.0)),
    };
    let nwk_parsed = nwk_read_str("((a:1,b:1)x:1,c:1)root:0;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let constraints = load_date_constraints(&dates, &graph, &names)?;
    Ok((graph, constraints))
  }

  #[test]
  fn test_pipeline_estimate_coalescent_tc_report_carries_the_skyline_solve() -> Result<(), Report> {
    let (graph, constraints) = dated_tree()?;
    let params = SkylineParams {
      n_points: 3,
      ..SkylineParams::default()
    };

    let node_times = TimetreeState::seed_from_values(&graph, &constraints).coalescent_node_times()?;
    let solve = optimize_skyline(&graph, &params, &node_times)?;
    let timescale = estimate_coalescent_tc(CoalescentMode::Skyline, &graph, &params, &node_times)?
      .expect("skyline mode yields a coalescent timescale");
    let report = timescale
      .report
      .expect("an inferred skyline carries a per-segment report");

    pretty_assert_array_eq!(solve.segment_boundaries, report.segment_boundaries);
    pretty_assert_array_eq!(solve.tc_schedule.values().clone(), timescale.schedule.values().clone());
    assert_eq!(Some(solve.log_likelihood.value()), report.log_likelihood);
    let band = report.band.expect("an inferred skyline carries a confidence band");
    pretty_assert_array_eq!(solve.tc_lower_bounds, band.lower);
    pretty_assert_array_eq!(solve.tc_upper_bounds, band.upper);

    Ok(())
  }
}
