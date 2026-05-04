use crate::cli::rtt_chart::{
  print_clock_regression_chart, write_clock_regression_chart_png, write_clock_regression_chart_svg,
};
use crate::commands::clock::args::{BranchSplitArgs, TreetimeClockArgs};
use crate::commands::clock::assign_dates::assign_dates;
use crate::commands::clock::clock_filter::clock_filter_inplace;
use crate::commands::clock::clock_graph::GraphClock;
use crate::commands::clock::clock_model::ClockModel;
use crate::commands::clock::clock_output::write_clock_model;
use crate::commands::clock::clock_regression::{ClockParams, estimate_clock_model_with_reroot_policy};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::clock::reroot::RerootParams;
use crate::commands::clock::rtt::{gather_clock_regression_results, write_clock_regression_result_csv};
use crate::make_report;
use eyre::{Report, WrapErr};
use log::info;
use treetime_io::dates_csv::read_dates;
use treetime_io::graphviz::graphviz_write_file;
use treetime_io::json::{JsonPretty, json_write_file};
use treetime_io::nwk::{NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime_utils::io::console::is_tty;

pub fn run_clock(clock_args: &TreetimeClockArgs) -> Result<(), Report> {
  let TreetimeClockArgs {
    aln,
    tree,
    vcf_reference,
    dates,
    name_column,
    date_column,
    sequence_length,
    gtr,
    gtr_params,
    branch_length_mode,
    method_anc,
    clock_filter,
    reroot,
    keep_root,
    prune_short,
    tip_slack,
    covariation,
    allow_negative_rate,
    plot_rtt,
    outdir,
    seed,
    branch_split,
    clock_regression,
  } = clock_args;

  let mut graph: GraphClock = if let Some(tree) = tree {
    nwk_read_file(tree)
  } else {
    unimplemented!("Tree inference is not implemented")
  }?;

  {
    let dates = read_dates(dates, name_column, date_column).wrap_err("When reading dates")?;
    assign_dates(&graph, &dates)?;
  }

  // Split workflow into a separate blocks depending whether covariation is used or not
  let (clock_model, new_outliers) = if *covariation {
    let seq_len = sequence_length
      .ok_or_else(|| make_report!("--sequence-length is required when --covariation is enabled"))? as f64;
    let tip_slack = tip_slack.unwrap_or(3.0);
    let overdispersion = 2.0; // TODO: empirical value for now, need to think of a better parameter than `tip_slack`
    let options = ClockParams {
      variance_factor: overdispersion / seq_len,
      variance_offset: 0.0,
      variance_offset_leaf: tip_slack * tip_slack / seq_len / seq_len,
    };

    estimate_clock_model_with_prefilter(
      &mut graph,
      &options,
      *keep_root,
      branch_split,
      *clock_filter,
      *allow_negative_rate,
    )?
  } else {
    estimate_clock_model_with_prefilter(
      &mut graph,
      &clock_regression.clock_params,
      *keep_root,
      branch_split,
      *clock_filter,
      *allow_negative_rate,
    )?
  };

  if let Some(delta) = new_outliers {
    info!("Clock filter changed outlier status for {delta} leaf nodes");
  }

  nwk_write_file(outdir.join("rerooted.nwk"), &graph, &NwkWriteOptions::default())?;
  json_write_file(outdir.join("graph_output.json"), &graph, JsonPretty(true))?;
  graphviz_write_file(outdir.join("graph_output.dot"), &graph)?;

  write_clock_model(&clock_model, &outdir.join("clock_model"))?;

  let results = gather_clock_regression_results(&graph, &clock_model);

  write_clock_regression_result_csv(&results, outdir.join("clock.csv"), b',')?;
  write_clock_regression_chart_svg(&results, &clock_model, outdir.join("clock.svg"))?;
  write_clock_regression_chart_png(&results, &clock_model, outdir.join("clock.png"))?;

  if is_tty() {
    print_clock_regression_chart(&results, &clock_model)?;
  }

  Ok(())
}

fn estimate_clock_model_with_prefilter(
  graph: &mut GraphClock,
  options: &ClockParams,
  keep_root: bool,
  branch_split: &BranchSplitArgs,
  clock_filter_threshold: f64,
  allow_negative_rate: bool,
) -> Result<(ClockModel, Option<i32>), Report> {
  let delta = (clock_filter_threshold > 0.0)
    .then(|| -> Result<i32, Report> {
      let params = BranchPointOptimizationParams::from(branch_split);
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
        &params,
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

  let params = BranchPointOptimizationParams::from(branch_split);
  let reroot_params = RerootParams {
    force_positive_rate: !allow_negative_rate,
    ..RerootParams::default()
  };
  let result = estimate_clock_model_with_reroot_policy(graph, options, None, keep_root, &params, &reroot_params, None)
    .wrap_err_with(|| {
      if delta.is_some() {
        "Clock model estimation failed after outlier filtering. \
         The pre-filter step removed outliers but the clock rate remains negative at all root positions."
          .to_owned()
      } else {
        "Clock model estimation failed".to_owned()
      }
    })?;
  Ok((result.clock_model, delta))
}
