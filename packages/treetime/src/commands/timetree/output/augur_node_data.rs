use crate::clock::clock_model::ClockModel;
use crate::partition::timetree::GraphTimetree;
use crate::partition::traits::BranchTopology;
use crate::timetree::confidence::NodeConfidenceInterval;
use eyre::Report;
use std::collections::BTreeMap;
use std::path::Path;
use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
use treetime_graph::node::{GraphNodeKey, Named};
use treetime_io::dates_csv::{DateConstraint, DatesMap};
use treetime_utils::datetime::year_fraction::year_fraction_to_datestring;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use treetime_utils::make_internal_report;
use util_augur_node_data_json::{
  AugurNodeDataJsonClock, AugurNodeDataJsonGeneratedBy, AugurNodeDataJsonRefine, AugurNodeDataJsonRefineMeta,
  AugurNodeDataJsonRefineNode,
};

/// Write augur-compatible node data JSON for the `timetree` command.
///
/// Produces the structure consumed by `augur export v2 --node-data`, equivalent
/// to `augur refine` output: top-level `clock` regression parameters, `alignment`
/// and `input_tree` paths, and per-node dates, branch lengths, and confidence
/// intervals.
///
/// Augur builds this JSON itself from TreeTime tree-node attributes (it never
/// calls the TreeTime CLI), so the contract is augur's `refine.py` traced against
/// TreeTime's final node-attribute state, not a v0 CLI oracle. Branch-field
/// semantics, traced through `TreeTime.run()`:
///
/// - `mutation_length` = TreeTime `node.mutation_length`: the ML-optimized branch
///   length in substitutions/site (`treeanc.py` sets `node.mutation_length =
///   node.branch_length` during branch-length optimization). This drives
///   `node_attrs.div` in `augur export v2` (cumulative `mutation_length`), so it
///   must be the divergence length. Source: `EdgeTimetree.base.branch_length`.
/// - `clock_length` = TreeTime `node.clock_length = up.time_before_present -
///   time_before_present` (`clock_tree.py:924`): the time-tree branch duration in
///   years. Equals `child.numdate - parent.numdate`. Source: difference of
///   `NodeTimetree.time`.
/// - `branch_length` = `clock_length`: `make_time_tree` overwrites
///   `node.branch_length = node.clock_length` (`clock_tree.py:925`) and augur
///   emits it unchanged under the default `--divergence-units=mutations-per-site`.
///   So `branch_length` carries the time value, identical to `clock_length`.
///
/// `dates` carries the parsed metadata date constraints used for `raw_date` (tips)
/// and `date_inferred`; the inferred `date` string derives from each node's
/// `numdate`.
///
/// When `mutation_counts` is `Some`, `mutation_length` is set to the per-edge
/// mutation count instead of the ML branch length (subs/site). `branch_length`
/// and `clock_length` remain time-valued (years) regardless.
pub fn build_augur_node_data_json(
  graph: &GraphTimetree,
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
    let node_guard = node.read_arc();
    let node_key = node_guard.key();
    let is_leaf = node_guard.is_leaf();
    let payload = node_guard.payload().read_arc();
    let node_name = payload
      .name()
      .map_or_else(|| format!("node_{}", node_key.as_usize()), |n| n.as_ref().to_owned());
    let numdate = payload.time;

    // Per-branch fields live on the parent edge. The root has no incoming branch:
    // augur emits placeholder values there, but `export v2` sets the root divergence
    // to 0 regardless and never consumes the root's branch fields, so we omit them.
    let (branch_length, clock_length, mutation_length) = match graph.node_parent(node_key)? {
      Some((parent_key, edge_key)) => {
        let mutation_length = if let Some(counts) = mutation_counts {
          Some(counts.get(&edge_key).copied().unwrap_or_default() as f64)
        } else {
          let edge = graph
            .get_edge(edge_key)
            .ok_or_else(|| make_internal_report!("Timetree node data: missing edge {edge_key:?}"))?;
          edge.read_arc().payload().read_arc().branch_length()
        };
        let clock_length = parent_time(graph, parent_key)?
          .zip(numdate)
          .map(|(parent, child)| child - parent);
        (clock_length.unwrap_or(0.0), clock_length, mutation_length)
      },
      None => (0.0, None, None),
    };

    let constraint: Option<&DateConstraint> = dates.and_then(|dates| dates.get(&node_name)).and_then(Option::as_ref);

    // date_inferred mirrors augur `not isinstance(node.raw_date_constraint, float)`:
    // a node carrying an exact (point) date is not inferred; range, uncertain, or
    // absent constraints are inferred.
    let date_inferred = !constraint.is_some_and(DateConstraint::is_exact);

    // raw_date is the original metadata string, set on tips only (augur sets it
    // for terminals from the metadata table).
    let raw_date = if is_leaf {
      constraint.map(|constraint| constraint.raw.clone())
    } else {
      None
    };

    let date = numdate.map(year_fraction_to_datestring);
    let num_date_confidence = ci_map.as_ref().and_then(|ci_map| ci_map.get(&node_key).copied());

    let confidence = payload.base.confidence;

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
  graph: &GraphTimetree,
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

/// Build the top-level `clock` object from the root-to-tip regression.
///
/// Mirrors augur `refine.py`: `rtt_Tmrca = -intercept / rate`, `cov` is the 2x2
/// regression covariance over `[rate, intercept]`, and `rate_std = sqrt(cov[0, 0])`
/// (the rate variance). `cov` and `rate_std` are present only when the clock rate
/// was estimated (absent for a fixed rate).
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

/// Inferred numeric date (`numdate`) of a node by key, used for the parent endpoint
/// of `clock_length = child.numdate - parent.numdate`.
fn parent_time(graph: &GraphTimetree, parent_key: GraphNodeKey) -> Result<Option<f64>, Report> {
  let parent = graph
    .get_node(parent_key)
    .ok_or_else(|| make_internal_report!("Timetree node data: missing parent node {parent_key:?}"))?;
  let time = parent.read_arc().payload().read_arc().time;
  Ok(time)
}

/// Confidence-interval lookup keyed by `GraphNodeKey` for stable lookup
/// independent of display naming, matching the auspice writer.
fn build_ci_map(intervals: &[NodeConfidenceInterval]) -> BTreeMap<GraphNodeKey, [f64; 2]> {
  intervals.iter().map(|ci| (ci.key, [ci.lower, ci.upper])).collect()
}
