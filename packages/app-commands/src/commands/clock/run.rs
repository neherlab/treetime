use crate::commands::clock::args::{BranchSplitArgs, OptimizationMethodCli, TreetimeClockArgs};
use crate::commands::shared::dates_input::read_input_dates;
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::shared::resolve_outputs::ResolveOutputs;
use crate::commands::shared::tree_input::read_input_tree;
use crate::rtt_chart::{write_clock_regression_chart_png, write_clock_regression_chart_svg};
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeDates};
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use app_output::table_output::table_write_file;
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use treetime::cancel::Cancel;
use treetime::clock::clock_model::ClockModel;
use treetime::clock::clock_regression::ClockVarianceParams;
use treetime::clock::clock_state::ClockInputs;
use treetime::clock::find_best_root::params::BranchPointOptimizationParams;
use treetime::clock::pipeline::{self, ClockInput, ClockParams};
use treetime::clock::rtt::ClockRegressionResult;
use treetime::make_report;
use treetime::progress::{LogSink, StageSink};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::io::json::{JsonPretty, json_write_file};

#[expect(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn run_clock(
  clock_args: &TreetimeClockArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<ClockResult, Report> {
  cancel.check()?;
  stages.report("Reading input", 0.0, "");

  let nwk_parsed = read_input_tree(&clock_args.tree, log)?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let input_order = leaf_order(&graph);

  let dates = read_input_dates(
    clock_args.metadata(),
    &clock_args.metadata_id,
    clock_args.date_column.date_column.as_deref(),
    &graph,
    &names,
    log,
  )?;

  let resolved = clock_args.resolve_outputs()?;

  let clock_params = if clock_args.covariation {
    let seq_len = clock_args
      .sequence_length
      .ok_or_else(|| make_report!("--sequence-length is required when --covariation is enabled"))?
      as f64;
    let tip_slack = clock_args.tip_slack.unwrap_or(3.0);
    let overdispersion = 2.0;
    ClockVarianceParams {
      variance_factor: overdispersion / seq_len,
      variance_offset: 0.0,
      variance_offset_leaf: tip_slack * tip_slack / seq_len / seq_len,
    }
  } else {
    clock_args.clock_regression.clock_params.clone().into()
  };

  let params = ClockParams {
    clock_params,
    clock_filter: clock_args.clock_filter,
    keep_root: clock_args.keep_root,
    allow_negative_rate: clock_args.allow_negative_rate,
    branch_params: branch_split_to_params(&clock_args.branch_split),
    reroot_spec: clock_args.reroot.spec(&graph, &names)?,
  };

  let input = ClockInput {
    graph,
    dates,
    branch_lengths,
  };

  let output = pipeline::run(&params, input, &names, cancel, stages, log).map_err(|err| err.into_report())?;
  let pipeline::ClockOutput {
    mut graph,
    inputs,
    divergences,
    outliers,
    clock_model,
    regression_results,
    names,
    branch_lengths,
  } = output;
  let topology_order = clock_args
    .topology_order
    .resolve_topology_order(&graph, &names, Some(input_order))?;
  topology_order.apply(&mut graph, &names, &branch_lengths)?;
  stages.report("Writing output", 0.8, "");

  write_clock_run_outputs(&resolved, &clock_model, &regression_results, log)?;

  let dates = ClockDates::new(&graph, &inputs, &outliers);
  write_clock_trees(&graph, &names, &branch_lengths, &divergences, &dates, &resolved, log)?;

  stages.report("Done", 1.0, "");
  Ok(ClockResult {
    clock_model,
    regression_results,
  })
}

pub struct ClockResult {
  pub clock_model: ClockModel,
  pub regression_results: Vec<ClockRegressionResult>,
}

fn write_clock_run_outputs(
  resolved: &ResolvedOutputs,
  clock_model: &ClockModel,
  regression_results: &[ClockRegressionResult],
  log: &dyn LogSink,
) -> Result<(), Report> {
  if let Some(path) = resolved.path(OutputSelection::ClockModel) {
    json_write_file(path, clock_model, JsonPretty(true))?;
  }
  if let Some(path) = resolved.path(OutputSelection::ClockCsv) {
    table_write_file(OutputSelection::ClockCsv, path, regression_results)?;
  }
  if let Some(path) = resolved.path(OutputSelection::ClockChartSvg) {
    write_clock_regression_chart_svg(regression_results, clock_model, path)?;
  }
  if let Some(path) = resolved.path(OutputSelection::ClockChartPng) {
    write_clock_regression_chart_png(regression_results, clock_model, path, log)?;
  }
  Ok(())
}

fn write_clock_trees(
  graph: &Graph,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  divergences: &BTreeMap<GraphNodeKey, f64>,
  dates: &ClockDates,
  resolved: &ResolvedOutputs,
  log: &dyn LogSink,
) -> Result<(), Report> {
  let annotated = AnnotatedGraph {
    graph,
    names,
    divergence_branch_lengths: branch_lengths,
    time_branch_lengths: None,
    divergence: Divergence::Values(divergences),
    sequences: None,
    dates: Some(TreeDates {
      num_date: &dates.num_date,
      confidence: None,
      excluded: &dates.excluded,
      input_dates: None,
    }),
    traits: None,
  };
  write_graph_outputs(&annotated, &resolved.tree_outputs)?;
  if let Some(tree) = tree_view_for_outputs(&annotated, resolved)? {
    write_tree_outputs(&tree, &resolved.tree_outputs, CommandKind::Clock, log)?;
  }
  Ok(())
}

struct ClockDates {
  num_date: BTreeMap<GraphNodeKey, Option<f64>>,
  excluded: BTreeSet<GraphNodeKey>,
}

impl ClockDates {
  fn new(graph: &Graph, inputs: &ClockInputs, outliers: &BTreeSet<GraphNodeKey>) -> Self {
    let num_date = graph
      .get_nodes()
      .map(|node| (node.key(), inputs.node(node.key()).time))
      .collect();
    let excluded = graph
      .get_nodes()
      .map(|node| node.key())
      .filter(|key| outliers.contains(key) || inputs.node(*key).bad_branch)
      .collect();
    Self { num_date, excluded }
  }
}

fn branch_split_to_params(args: &BranchSplitArgs) -> BranchPointOptimizationParams {
  match args.method {
    OptimizationMethodCli::Grid => BranchPointOptimizationParams::grid_with(args.grid_params.clone().into()),
    OptimizationMethodCli::Brent => BranchPointOptimizationParams::brent_with(args.brent_params.clone().into()),
    OptimizationMethodCli::GoldenSection => {
      BranchPointOptimizationParams::golden_section_with(args.golden_params.clone().into())
    },
  }
}
