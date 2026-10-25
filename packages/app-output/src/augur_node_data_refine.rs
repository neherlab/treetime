use crate::annotated_graph::{AnnotatedTreeView, TreeDates};
use eyre::Report;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use treetime::clock::clock_model::ClockModel;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::key_map::EdgeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::DateConstraint;
use treetime_utils::cast::ToF64Exact;
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
      let name = graph.names[key].to_owned();
      let lengths = refine_lengths(tree, run, key);
      let dates = graph
        .dates
        .as_ref()
        .map(|dates| refine_dates(dates, key, tree.tree().children(key).is_empty()))
        .unwrap_or_default();
      let node = AugurNodeDataJsonRefineNode {
        branch_length: lengths.branch_length,
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
      alignment: (!run.alignments.is_empty()).then(|| run.alignments.iter().map(|path| path.display()).join(",")),
      input_tree: run.input_tree.map(|path| path.display().to_string()),
      clock: run.clock_model.map(build_clock),
      other: BTreeMap::new(),
    },
    nodes,
  }
}

pub fn total_mutation_counts<'a>(counts: impl IntoIterator<Item = &'a EdgeMap<usize>>) -> EdgeMap<usize> {
  counts
    .into_iter()
    .flatten()
    .fold(EdgeMap::new(), |mut total, (edge, &count)| {
      *total.entry(edge).or_default() += count;
      total
    })
}

fn refine_lengths(tree: &AnnotatedTreeView<'_>, run: &RefineRun<'_>, key: GraphNodeKey) -> RefineLengths {
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
        mutation_length: mutation_length(tree, run, edge_key),
      }
    },
    (None, edge_key) => RefineLengths {
      branch_length: edge_key
        .and_then(|edge_key| mutation_length(tree, run, edge_key))
        .unwrap_or(0.0),
      clock_length: None,
      mutation_length: None,
    },
  }
}

fn mutation_length(tree: &AnnotatedTreeView<'_>, run: &RefineRun<'_>, edge_key: GraphEdgeKey) -> Option<f64> {
  match run.mutation_counts {
    Some(counts) => Some(counts.get(edge_key).copied().unwrap_or_default().to_f64_exact()),
    None => tree.graph().divergence_branch_lengths[&edge_key],
  }
}

pub struct RefineRun<'a> {
  pub alignments: &'a [PathBuf],
  pub mutation_counts: Option<&'a EdgeMap<usize>>,
  pub input_tree: Option<&'a Path>,
  pub clock_model: Option<&'a ClockModel>,
}

fn refine_dates(dates: &TreeDates<'_>, key: GraphNodeKey, is_leaf: bool) -> RefineDates {
  let constraint = dates.input_dates.and_then(|input_dates| input_dates.get(key));
  let numdate = dates.num_date[&key];
  RefineDates {
    numdate,
    date: numdate.map(year_fraction_to_datestring),
    num_date_confidence: dates.confidence.and_then(|confidence| confidence.get(key).copied()),
    date_inferred: Some(!constraint.is_some_and(DateConstraint::is_exact)),
    raw_date: constraint.filter(|_| is_leaf).map(|constraint| constraint.raw.clone()),
  }
}

fn build_clock(clock_model: &ClockModel) -> AugurNodeDataJsonClock {
  let rate = clock_model.clock_rate();
  let intercept = clock_model.intercept();
  AugurNodeDataJsonClock {
    rate,
    intercept,
    rtt_tmrca: -intercept / rate,
    cov: None,
    rate_std: clock_model.rate_std(),
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
