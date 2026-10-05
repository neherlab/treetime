use crate::annotated_graph::{AnnotatedTreeView, TreeDates};
use eyre::Report;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::clock::clock_model::ClockModel;
use treetime_graph::assign_node_names::node_name_or_key;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::DateConstraint;
use treetime_utils::datetime::year_fraction::year_fraction_to_datestring;
use treetime_utils::io::json::{JsonPretty, json_write_file};
use util_augur_node_data_json::{
  AugurNodeDataJsonClock, AugurNodeDataJsonGeneratedBy, AugurNodeDataJsonRefine, AugurNodeDataJsonRefineMeta,
  AugurNodeDataJsonRefineNode,
};

pub fn write_augur_node_data_refine(
  tree: &AnnotatedTreeView<'_>,
  run: &RefineRun<'_>,
  path: &Path,
) -> Result<(), Report> {
  json_write_file(path, &build_augur_node_data_refine(tree, run), JsonPretty(true))
}

pub fn build_augur_node_data_refine(tree: &AnnotatedTreeView<'_>, run: &RefineRun<'_>) -> AugurNodeDataJsonRefine {
  let graph = tree.graph();
  let nodes = graph
    .graph
    .get_nodes()
    .map(|node| {
      let key = node.key();
      let name = node_name_or_key(key, graph.names[&key].as_deref());
      let lengths = refine_lengths(tree, key);
      let dates = graph
        .dates
        .as_ref()
        .map(|dates| refine_dates(dates, key, &name, tree.tree().children(key).is_empty()))
        .unwrap_or_default();
      let node = AugurNodeDataJsonRefineNode {
        branch_length: lengths.branch_length,
        confidence: run
          .branch_support
          .and_then(|support| support.get(&key).copied().flatten()),
        numdate: dates.numdate,
        clock_length: lengths.clock_length,
        mutation_length: lengths.mutation_length,
        raw_date: dates.raw_date,
        date: dates.date,
        date_inferred: dates.date_inferred,
        num_date_confidence: dates.num_date_confidence,
        other: BTreeMap::new(),
      };
      (name, node)
    })
    .collect();

  AugurNodeDataJsonRefine {
    generated_by: Some(AugurNodeDataJsonGeneratedBy {
      program: "treetime".to_owned(),
      version: env!("CARGO_PKG_VERSION").to_owned(),
    }),
    metadata: AugurNodeDataJsonRefineMeta {
      alignment: run.alignment.map(|path| path.display().to_string()),
      input_tree: run.input_tree.map(|path| path.display().to_string()),
      clock: run.clock_model.map(build_clock),
      other: BTreeMap::new(),
    },
    nodes,
  }
}

pub struct RefineRun<'a> {
  pub alignment: Option<&'a Path>,
  pub input_tree: Option<&'a Path>,
  pub clock_model: Option<&'a ClockModel>,
  pub branch_support: Option<&'a BTreeMap<GraphNodeKey, Option<f64>>>,
}

fn refine_lengths(tree: &AnnotatedTreeView<'_>, key: GraphNodeKey) -> RefineLengths {
  let graph = tree.graph();
  let parent_edge = tree.tree().parent(key).map(|(_, edge_key)| edge_key);
  match (graph.time_branch_lengths, parent_edge) {
    (Some(_), None) => RefineLengths {
      branch_length: 0.0,
      clock_length: Some(0.0),
      mutation_length: Some(0.0),
    },
    (Some(time_lengths), Some(edge_key)) => {
      let clock_length = time_lengths[&edge_key];
      RefineLengths {
        branch_length: clock_length.unwrap_or(0.0),
        clock_length,
        mutation_length: mutation_length(tree, edge_key),
      }
    },
    (None, edge_key) => RefineLengths {
      branch_length: edge_key
        .and_then(|edge_key| mutation_length(tree, edge_key))
        .unwrap_or(0.0),
      clock_length: None,
      mutation_length: None,
    },
  }
}

#[expect(
  clippy::as_conversions,
  reason = "a mutation count is far below 2^53, so the conversion to f64 is exact"
)]
fn mutation_length(tree: &AnnotatedTreeView<'_>, edge_key: GraphEdgeKey) -> Option<f64> {
  let graph = tree.graph();
  match graph.sequences.as_ref().and_then(|sequences| sequences.mutation_counts) {
    Some(counts) => Some(counts.get(&edge_key).copied().unwrap_or_default() as f64),
    None => graph.divergence_branch_lengths[&edge_key],
  }
}

fn refine_dates(dates: &TreeDates<'_>, key: GraphNodeKey, name: &str, is_leaf: bool) -> RefineDates {
  let constraint = dates
    .input_dates
    .and_then(|input_dates| input_dates.get(name))
    .and_then(Option::as_ref);
  let numdate = dates.num_date[&key];
  RefineDates {
    numdate,
    date: numdate.map(year_fraction_to_datestring),
    num_date_confidence: dates.confidence.and_then(|confidence| confidence.get(&key).copied()),
    date_inferred: Some(!constraint.is_some_and(DateConstraint::is_exact)),
    raw_date: constraint.filter(|_| is_leaf).map(|constraint| constraint.raw.clone()),
  }
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

struct RefineLengths {
  branch_length: f64,
  clock_length: Option<f64>,
  mutation_length: Option<f64>,
}

#[derive(Default)]
struct RefineDates {
  numdate: Option<f64>,
  date: Option<String>,
  num_date_confidence: Option<[f64; 2]>,
  date_inferred: Option<bool>,
  raw_date: Option<String>,
}
