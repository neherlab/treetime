#![allow(
  clippy::as_conversions,
  reason = "application layer: counts and indices to f64, integer division for averaging, graph node access and CLI/config setup invariants, and default variant matches"
)]

//! N-API `clock` request shape and orchestration.

use crate::commands::support::{default_output_plan, default_topology_order, reroot_spec};
use app_output::clock_result::{ClockNodeOut, EdgeOut};
use app_output::clock_tree_output::write_clock_tree_outputs;
use app_output::output_plan::{CommandKind, OutputSelection};
use app_output::rtt::write_clock_regression_result_csv;
use eyre::{Report, WrapErr};
use serde::Deserialize;
use smart_default::SmartDefault;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::ancestral::params::MethodAncestral;
use treetime::cancel::Cancel;
use treetime::clock::clock_model::ClockModel;
use treetime::clock::clock_output::write_clock_model;
use treetime::clock::clock_regression::ClockVarianceParams;
use treetime::clock::clock_state::{ClockInputs, ClockState};
use treetime::clock::find_best_root::params::{BranchPointOptimizationParams, RerootMethod};
use treetime::clock::pipeline::{self, ClockInput, ClockParams};
use treetime::clock::rtt::ClockRegressionResult;
use treetime::gtr::get_gtr::GtrModelName;
use treetime::make_error;
use treetime::make_report;
use treetime::optimize::params::BranchLengthMode;
use treetime::progress::ProgressSink;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::csv::{default_metadata_delimiters, default_name_candidates};
use treetime_io::dates_csv::read_dates;
use treetime_io::nwk::{CommentProviders, nwk_read_file};

/// Clock estimation request (openapi subset).
#[derive(Debug, SmartDefault, Deserialize)]
#[serde(default)]
pub struct ClockArgs {
  pub aln: Vec<String>,
  pub tree: Option<String>,
  pub vcf_reference: Option<String>,
  pub dates: String,
  pub name_column: Option<String>,
  pub date_column: Option<String>,
  pub sequence_length: Option<usize>,
  #[default(GtrModelName::default())]
  pub gtr: GtrModelName,
  pub gtr_params: Vec<String>,
  #[default(BranchLengthMode::default())]
  pub branch_length_mode: BranchLengthMode,
  #[default(MethodAncestral::default())]
  pub method_anc: MethodAncestral,
  #[default = 3.0]
  pub clock_filter: f64,
  pub reroot: Option<RerootMethod>,
  pub reroot_tips: Vec<String>,
  pub keep_root: bool,
  pub prune_short: bool,
  pub tip_slack: Option<f64>,
  pub covariation: bool,
  pub allow_negative_rate: bool,
  pub outdir: String,
  pub seed: Option<u64>,
}

/// Clock estimation result. Every field is a value the caller may read; the serialized response body
/// carries none of them (they are internal to the run).
#[derive(serde::Serialize)]
pub struct ClockResult {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub nodes: BTreeMap<GraphNodeKey, ClockNodeOut>,
  #[serde(skip)]
  pub edges: BTreeMap<GraphEdgeKey, EdgeOut>,
  #[serde(skip)]
  pub clock_model: ClockModel,
  #[serde(skip)]
  pub regression_results: Vec<ClockRegressionResult>,
}

pub fn run_clock(args: &ClockArgs, cancel: &dyn Cancel, progress: &dyn ProgressSink) -> Result<ClockResult, Report> {
  cancel.check()?;
  progress.report("Reading input", 0.0, "");

  let nwk_parsed = if let Some(tree) = &args.tree {
    nwk_read_file(Path::new(tree))
  } else {
    return make_error!("Tree inference is not implemented. Provide a tree file with --tree");
  }?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let branch_lengths = nwk_parsed.branch_lengths;

  let id_columns = args
    .name_column
    .clone()
    .map_or_else(default_name_candidates, |col| vec![col]);
  let dates = read_dates(
    Path::new(&args.dates),
    &default_metadata_delimiters(),
    &id_columns,
    &None,
    &args.date_column,
  )
  .wrap_err("When reading dates")?;

  let resolved = default_output_plan(CommandKind::Clock, Path::new(&args.outdir))?;

  let clock_params = if args.covariation {
    let seq_len = args
      .sequence_length
      .ok_or_else(|| make_report!("--sequence-length is required when --covariation is enabled"))?
      as f64;
    let tip_slack = args.tip_slack.unwrap_or(3.0);
    let overdispersion = 2.0;
    ClockVarianceParams {
      variance_factor: overdispersion / seq_len,
      variance_offset: 0.0,
      variance_offset_leaf: tip_slack * tip_slack / seq_len / seq_len,
    }
  } else {
    ClockVarianceParams::default()
  };

  let params = ClockParams {
    clock_params,
    clock_filter: args.clock_filter,
    keep_root: args.keep_root,
    allow_negative_rate: args.allow_negative_rate,
    branch_params: BranchPointOptimizationParams::default(),
    reroot_spec: reroot_spec(args.reroot, &args.reroot_tips),
  };

  let input = ClockInput {
    graph,
    dates,
    branch_lengths,
  };

  let output = pipeline::run(&params, input, &names, cancel, progress).map_err(|err| err.into_report())?;
  let pipeline::ClockOutput {
    mut graph,
    inputs,
    state,
    clock_model,
    regression_results,
    names,
    branch_lengths,
  } = output;
  default_topology_order().apply(&mut graph, &names, &branch_lengths)?;
  progress.report("Writing output", 0.8, "");

  let (nodes, edges) = gather_clock_outputs(&graph, &inputs, &state, &names, &branch_lengths);

  if !resolved.tree_outputs.is_empty() {
    write_clock_tree_outputs(
      &graph,
      &nodes,
      &branch_lengths,
      &resolved.tree_outputs,
      &CommentProviders::new(),
    )?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockModel) {
    write_clock_model(&clock_model, path)?;
  }

  if let Some(path) = resolved.non_tree_outputs.get(&OutputSelection::ClockCsv) {
    write_clock_regression_result_csv(&regression_results, path, b',')?;
  }

  progress.report("Done", 1.0, "");
  Ok(ClockResult {
    graph,
    nodes,
    edges,
    clock_model,
    regression_results,
  })
}

fn gather_clock_outputs(
  graph: &Graph,
  inputs: &ClockInputs,
  state: &ClockState,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
) -> (BTreeMap<GraphNodeKey, ClockNodeOut>, BTreeMap<GraphEdgeKey, EdgeOut>) {
  let nodes = graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      let name = names[&key].clone();
      let node_state = state.node(key);
      let node_input = inputs.node(key);
      let out = ClockNodeOut {
        name,
        div: node_state.div,
        time: node_input.time,
        is_outlier: node_state.is_outlier,
        bad_branch: node_input.bad_branch,
      };
      (key, out)
    })
    .collect();

  let edges = graph
    .get_edges()
    .map(|edge| {
      let key = edge.key();
      (
        key,
        EdgeOut {
          branch_length: branch_lengths[&key],
        },
      )
    })
    .collect();

  (nodes, edges)
}
