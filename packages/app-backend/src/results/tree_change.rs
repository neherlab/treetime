use crate::results::tree::ResultTree;
use deser::Serialize;
use eyre::Report;
use schemars::JsonSchema;
use std::path::Path;
use treetime_io::nwk::{NewickDialect, NwkParse};
use treetime_io::tree::tree_read_file;
use treetime_schema::skip_serializing_optionals;
use treetime_utils::cast::ToF64Exact;

pub fn tree_change(before: Option<&TreeStatistics>, after: Option<&ResultTree>) -> TreeChange {
  let mutations = after.map_or(0, |tree| tree.nodes.iter().map(|node| node.mutations.len()).sum());
  let after = after.map(result_tree_statistics);
  TreeChange {
    rows: statistic_rows(before, after.as_ref()),
    mutations,
  }
}

/// The tree before and after an `optimize` or `prune` run.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct TreeChange {
  /// Number of nucleotide mutations on all branches of the tree the run wrote.
  pub mutations: usize,
  /// The statistics of both trees by name, in the order to show them.
  pub rows: Vec<TreeStatisticRow>,
}

/// One statistic of the input tree and the result tree.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct TreeStatisticRow {
  /// Name of the statistic.
  pub label: String,
  /// Whether the statistic counts nodes; otherwise it is a length in substitutions per site.
  pub count: bool,
  /// Value for the input tree; absent when it cannot be read, and the problems of the results name the cause.
  pub before: Option<f64>,
  /// Value for the tree the run wrote; absent when the run wrote no Auspice tree.
  pub after: Option<f64>,
}

/// Size and shape of a tree.
#[derive(Clone, Debug, PartialEq)]
pub struct TreeStatistics {
  /// Number of samples.
  pub samples: usize,
  /// Number of internal nodes, the root included.
  pub internal_nodes: usize,
  /// Number of internal nodes with more than two children.
  pub polytomies: usize,
  /// Sum of the branch lengths, in substitutions per site; absent when the tree has no branch lengths.
  pub total_branch_length: Option<f64>,
}

pub fn input_tree_statistics(path: &Path, dialect: NewickDialect) -> Result<TreeStatistics, Report> {
  let NwkParse {
    graph, branch_lengths, ..
  } = tree_read_file(path, dialect)?;
  let lengths = branch_lengths.values().flatten().copied().collect::<Vec<_>>();
  Ok(TreeStatistics {
    samples: graph.num_leaves(),
    internal_nodes: graph.num_nodes() - graph.num_leaves(),
    polytomies: graph.get_nodes().filter(|node| node.outbound().len() > 2).count(),
    total_branch_length: (!lengths.is_empty()).then(|| lengths.iter().sum()),
  })
}

fn statistic_rows(before: Option<&TreeStatistics>, after: Option<&TreeStatistics>) -> Vec<TreeStatisticRow> {
  let row = |label: &str, count: bool, value: fn(&TreeStatistics) -> Option<f64>| TreeStatisticRow {
    label: label.to_owned(),
    count,
    before: before.and_then(value),
    after: after.and_then(value),
  };
  vec![
    row("Samples", true, |statistics| Some(statistics.samples.to_f64_exact())),
    row("Internal nodes", true, |statistics| {
      Some(statistics.internal_nodes.to_f64_exact())
    }),
    row("Polytomies", true, |statistics| {
      Some(statistics.polytomies.to_f64_exact())
    }),
    row("Total branch length", false, |statistics| {
      statistics.total_branch_length
    }),
  ]
}

fn result_tree_statistics(tree: &ResultTree) -> TreeStatistics {
  let samples = tree.tips().count();
  let lengths = tree
    .nodes
    .iter()
    .filter_map(|node| Some(node.div? - tree.nodes[node.parent?].div?))
    .collect::<Vec<_>>();
  TreeStatistics {
    samples,
    internal_nodes: tree.nodes.len() - samples,
    polytomies: tree.nodes.iter().filter(|node| node.children.len() > 2).count(),
    total_branch_length: (!lengths.is_empty()).then(|| lengths.iter().sum()),
  }
}
