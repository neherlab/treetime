use crate::results::tree::ResultTree;
use deser::Serialize;
use eyre::Report;
use schemars::JsonSchema;
use std::path::Path;
use treetime_io::nwk::{NewickDialect, NwkParse};
use treetime_io::tree::tree_read_file;
use treetime_schema::skip_serializing_optionals;

pub fn tree_change(before: Option<TreeStatistics>, after: Option<&ResultTree>) -> TreeChange {
  TreeChange {
    before,
    after: after.map(result_tree_statistics),
    mutations: after.map_or(0, |tree| tree.nodes.iter().map(|node| node.mutations.len()).sum()),
  }
}

/// The tree before and after an `optimize` or `prune` run.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
pub struct TreeChange {
  /// Statistics of the input tree; absent when it cannot be read, and the problems of the results name the cause.
  pub before: Option<TreeStatistics>,
  /// Statistics of the tree the run wrote; absent when the run wrote no Auspice tree.
  pub after: Option<TreeStatistics>,
  /// Number of nucleotide mutations on all branches of the tree the run wrote.
  pub mutations: usize,
}

/// Size and shape of a tree.
#[derive(Clone, Debug, PartialEq, JsonSchema, Serialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
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
