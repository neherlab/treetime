use crate::clock::assign_dates::assign_dates;
use crate::clock::clock_filter::clock_filter_inplace;
use crate::clock::clock_graph::GraphClock;
use crate::clock::clock_model::ClockModel;
use crate::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::clock::reroot::RerootParams;
use crate::clock::rtt::{ClockRegressionResult, gather_clock_regression_results};
use crate::progress::ProgressSink;
use eyre::{Report, WrapErr};
use log::info;
use serde::Serialize;
use treetime_io::dates_csv::DatesMap;

pub struct ClockPipelineParams {
  pub clock_params: ClockParams,
  pub clock_filter: f64,
  pub keep_root: bool,
  pub allow_negative_rate: bool,
  pub branch_params: BranchPointOptimizationParams,
}

pub struct ClockInput {
  pub graph: GraphClock,
  pub dates: DatesMap,
}

#[derive(Debug, Serialize)]
pub struct ClockOutput {
  #[serde(skip)]
  pub graph: GraphClock,
  pub clock_model: ClockModel,
  pub regression_results: Vec<ClockRegressionResult>,
}

pub fn run(
  params: &ClockPipelineParams,
  mut input: ClockInput,
  progress: &dyn ProgressSink,
) -> Result<ClockOutput, Report> {
  progress.check_cancelled()?;
  progress.report("Assigning dates", 0.1, "");
  assign_dates(&input.graph, &input.dates)?;

  progress.check_cancelled()?;
  progress.report("Clock regression", 0.3, "");
  let (clock_model, new_outliers) = estimate_clock_model_with_prefilter(
    &mut input.graph,
    &params.clock_params,
    params.keep_root,
    &params.branch_params,
    params.clock_filter,
    params.allow_negative_rate,
  )?;

  if let Some(delta) = new_outliers {
    info!("Clock filter changed outlier status for {delta} leaf nodes");
  }

  let regression_results = gather_clock_regression_results(&input.graph, &clock_model);

  progress.report("Done", 1.0, "");
  Ok(ClockOutput {
    graph: input.graph,
    clock_model,
    regression_results,
  })
}

fn estimate_clock_model_with_prefilter(
  graph: &mut GraphClock,
  options: &ClockParams,
  keep_root: bool,
  branch_params: &BranchPointOptimizationParams,
  clock_filter_threshold: f64,
  allow_negative_rate: bool,
) -> Result<(ClockModel, Option<i32>), Report> {
  let delta = (clock_filter_threshold > 0.0)
    .then(|| -> Result<i32, Report> {
      // Allow negative rates during pre-filter root finding. Some datasets (e.g. dengue/100)
      // have negative estimated rate at ALL root positions when outliers are included.
      // The pre-filter clock model only needs to be good enough for IQD-based outlier detection.
      let reroot_params = RerootParams {
        force_positive_rate: false,
        ..RerootParams::default()
      };
      let result = estimate_clock_model_with_reroot_policy(
        graph,
        &ClockParams::default(),
        None,
        keep_root,
        branch_params,
        &reroot_params,
        None,
      )?;
      let pre_clock_model = result.clock_model;
      if pre_clock_model.clock_rate() < 0.0 {
        // IQD-based filtering uses |deviation| > IQD * threshold, so the absolute-value
        // comparison is slope-sign-invariant: outliers are identified by distance from the
        // fitted line regardless of slope direction.
        log::warn!(
          "Pre-filter clock rate is negative ({:.6e}). Outlier detection proceeds with this model.",
          pre_clock_model.clock_rate()
        );
      }
      Ok(clock_filter_inplace(graph, &pre_clock_model, clock_filter_threshold).new_outliers)
    })
    .transpose()?;

  let reroot_params = RerootParams {
    force_positive_rate: !allow_negative_rate,
    ..RerootParams::default()
  };
  let result = estimate_clock_model_with_reroot_policy(graph, options, None, keep_root, branch_params, &reroot_params, None)
    .wrap_err_with(|| {
      if delta.is_some() {
        "Clock model estimation failed after outlier filtering. The pre-filter step removed outliers but the clock rate remains negative at all root positions.".to_owned()
      } else {
        "Clock model estimation failed".to_owned()
      }
    })?;
  Ok((result.clock_model, delta))
}
