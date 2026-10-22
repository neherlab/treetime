use crate::commands::clock::args::TreetimeClockArgs;
use crate::commands::shared::dates_input::read_input_dates;
use crate::commands::shared::leaf_order::leaf_order;
use crate::commands::shared::output_args::ResolveOutputs;
use crate::commands::shared::tree_input::read_input_tree;
#[cfg(not(target_family = "wasm"))]
use crate::rtt_chart::write_clock_regression_charts;
use app_output::annotated_graph::{AnnotatedGraph, Divergence, TreeDates};
#[cfg(target_family = "wasm")]
use app_output::output_plan::output_unavailable;
use app_output::output_plan::{CommandKind, OutputSelection, ResolvedOutputs};
use app_output::table_output::table_write_file;
use app_output::tree_output::{tree_view_for_outputs, write_graph_outputs, write_tree_outputs};
use eyre::Report;
#[cfg(not(target_family = "wasm"))]
use itertools::Itertools;
use std::collections::{BTreeMap, BTreeSet};
use treetime::cancel::Cancel;
use treetime::clock::clock_model::ClockModel;
use treetime::clock::clock_state::ClockInputs;
use treetime::clock::pipeline::{self, ClockInput, ClockParams};
use treetime::clock::rtt::ClockRegressionResult;
use treetime::progress::{LogSink, StageSink};
#[cfg(not(target_family = "wasm"))]
use treetime_charts::convert::ChartFormat;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_utils::io::json::{JsonPretty, json_write_file};

pub fn run_clock(
  clock_args: &TreetimeClockArgs,
  cancel: &dyn Cancel,
  stages: &dyn StageSink,
  log: &dyn LogSink,
) -> Result<ClockResult, Report> {
  cancel.check()?;
  stages.report("Reading input", 0.0, "");

  let nwk_parsed = read_input_tree(&clock_args.tree, clock_args.tree_dialect.dialect(), log)?;
  let names = nwk_parsed.names;
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;
  let input_order = leaf_order(&graph);

  let dates = read_input_dates(
    clock_args.metadata(),
    &clock_args.metadata_id,
    clock_args.date_column.date_column.as_deref(),
    &clock_args.date_column.date_format,
    &graph,
    &names,
    log,
  )?;

  let resolved = clock_args.resolve_outputs()?;

  let params = ClockParams {
    clock_filter: clock_args.clock_filter,
    keep_root: clock_args.keep_root,
    allow_negative_rate: clock_args.allow_negative_rate,
    reroot_spec: clock_args.reroot.spec(&graph, &names)?,
  };

  let input = ClockInput {
    graph,
    names,
    dates,
    branch_lengths,
  };

  let output = pipeline::run(&params, input, cancel, stages, log).map_err(|err| err.into_report())?;
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
  let topology_order = clock_args.topology_order.resolve_topology_order(
    &graph,
    &names,
    Some(input_order),
    clock_args.tree_dialect.dialect(),
  )?;
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
  write_clock_charts(resolved, clock_model, regression_results, log)
}

#[cfg(not(target_family = "wasm"))]
fn write_clock_charts(
  resolved: &ResolvedOutputs,
  clock_model: &ClockModel,
  regression_results: &[ClockRegressionResult],
  log: &dyn LogSink,
) -> Result<(), Report> {
  let charts = [
    (ChartFormat::Svg, OutputSelection::ClockChartSvg),
    (ChartFormat::Png, OutputSelection::ClockChartPng),
    (ChartFormat::Pdf, OutputSelection::ClockChartPdf),
  ]
  .into_iter()
  .filter_map(|(format, selection)| resolved.path(selection).map(|path| (format, path)))
  .collect_vec();
  write_clock_regression_charts(regression_results, clock_model, &charts)
}

#[cfg(target_family = "wasm")]
fn write_clock_charts(
  resolved: &ResolvedOutputs,
  clock_model: &ClockModel,
  regression_results: &[ClockRegressionResult],
  log: &dyn LogSink,
) -> Result<(), Report> {
  [
    OutputSelection::ClockChartSvg,
    OutputSelection::ClockChartPng,
    OutputSelection::ClockChartPdf,
  ]
  .into_iter()
  .filter_map(|selection| resolved.non_tree_outputs.get(&selection).map(|file| (selection, file)))
  .try_for_each(|(selection, file)| output_unavailable(selection, file, "the web app back end draws no charts", log))
}

fn write_clock_trees(
  graph: &Graph,
  names: &NodeNames,
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
