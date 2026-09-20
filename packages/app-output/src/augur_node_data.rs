use crate::timetree_result::{TimetreeEdgeOut, TimetreeNodeOut};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::clock::clock_model::ClockModel;
use treetime::timetree::confidence::NodeConfidenceInterval;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::{DateConstraint, DatesMap};
use treetime_utils::datetime::year_fraction::year_fraction_to_datestring;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use util_augur_node_data_json::{
  AugurNodeDataJsonClock, AugurNodeDataJsonGeneratedBy, AugurNodeDataJsonRefine, AugurNodeDataJsonRefineMeta,
  AugurNodeDataJsonRefineNode,
};

#[allow(
  clippy::as_conversions,
  reason = "count/index numeric cast is exact for the domain range"
)]
pub fn build_augur_node_data_json(
  graph: &Graph,
  outputs: &BTreeMap<GraphNodeKey, TimetreeNodeOut>,
  edges: &BTreeMap<GraphEdgeKey, TimetreeEdgeOut>,
  clock_model: &ClockModel,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  dates: Option<&DatesMap>,
  alignment: Option<&Path>,
  input_tree: Option<&Path>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
) -> Result<AugurNodeDataJsonRefine, Report> {
  let ci_map = confidence_intervals.map(build_ci_map);

  let mut nodes = BTreeMap::new();
  for node in graph.get_nodes() {
    let node_guard = node;
    let node_key = node_guard.key();
    let out = &outputs[&node_key];
    let is_leaf = node_guard.is_leaf();
    let node_name = out
      .name
      .as_deref()
      .map_or_else(|| format!("node_{}", node_key.as_usize()), str::to_owned);
    let numdate = out.time;

    let (branch_length, clock_length, mutation_length) = match graph.node_parent(node_key)? {
      Some((parent_key, edge_key)) => {
        let mutation_length = if let Some(counts) = mutation_counts {
          Some(counts.get(&edge_key).copied().unwrap_or_default() as f64)
        } else {
          edges[&edge_key].branch_length
        };
        let clock_length = outputs[&parent_key]
          .time
          .zip(numdate)
          .map(|(parent, child)| child - parent);
        (clock_length.unwrap_or(0.0), clock_length, mutation_length)
      },
      None => (0.0, Some(0.0), Some(0.0)),
    };

    let constraint: Option<&DateConstraint> = dates.and_then(|dates| dates.get(&node_name)).and_then(Option::as_ref);

    let date_inferred = !constraint.is_some_and(DateConstraint::is_exact);

    let raw_date = if is_leaf {
      constraint.map(|constraint| constraint.raw.clone())
    } else {
      None
    };

    let date = numdate.map(year_fraction_to_datestring);
    let num_date_confidence = ci_map.as_ref().and_then(|ci_map| ci_map.get(&node_key).copied());

    let confidence = out.confidence;

    nodes.insert(
      node_name,
      AugurNodeDataJsonRefineNode {
        branch_length,
        confidence,
        numdate,
        clock_length,
        mutation_length,
        raw_date,
        date,
        date_inferred: Some(date_inferred),
        num_date_confidence,
        other: BTreeMap::new(),
      },
    );
  }

  Ok(AugurNodeDataJsonRefine {
    generated_by: Some(AugurNodeDataJsonGeneratedBy {
      program: "treetime".to_owned(),
      version: env!("CARGO_PKG_VERSION").to_owned(),
    }),
    metadata: AugurNodeDataJsonRefineMeta {
      alignment: alignment.map(|path| path.display().to_string()),
      input_tree: input_tree.map(|path| path.display().to_string()),
      clock: Some(build_clock(clock_model)),
      other: BTreeMap::new(),
    },
    nodes,
  })
}

pub fn write_augur_node_data_json(
  graph: &Graph,
  outputs: &BTreeMap<GraphNodeKey, TimetreeNodeOut>,
  edges: &BTreeMap<GraphEdgeKey, TimetreeEdgeOut>,
  clock_model: &ClockModel,
  confidence_intervals: Option<&[NodeConfidenceInterval]>,
  dates: Option<&DatesMap>,
  alignment: Option<&Path>,
  input_tree: Option<&Path>,
  mutation_counts: Option<&BTreeMap<GraphEdgeKey, usize>>,
  path: &Path,
) -> Result<(), Report> {
  let data = build_augur_node_data_json(
    graph,
    outputs,
    edges,
    clock_model,
    confidence_intervals,
    dates,
    alignment,
    input_tree,
    mutation_counts,
  )?;
  json_write_file(path, &data, JsonPretty(true))?;
  Ok(())
}

fn build_clock(clock_model: &ClockModel) -> AugurNodeDataJsonClock {
  let rate = clock_model.clock_rate();
  let intercept = clock_model.intercept();
  let cov = clock_model
    .cov()
    .map(|cov| cov.outer_iter().map(|row| row.to_vec()).collect());
  let rate_std = clock_model.cov().map(|cov| cov[[0, 0]].sqrt());

  AugurNodeDataJsonClock {
    rate,
    intercept,
    rtt_tmrca: -intercept / rate,
    cov,
    rate_std,
    other: BTreeMap::new(),
  }
}

fn build_ci_map(intervals: &[NodeConfidenceInterval]) -> BTreeMap<GraphNodeKey, [f64; 2]> {
  intervals.iter().map(|ci| (ci.key, [ci.lower, ci.upper])).collect()
}
