use crate::command::{AppCommand, OutputFile};
use crate::job::JobId;
use crate::results::outputs::read_auspice;
use crate::results::tree::{DateInterval, ResultTree};
use crate::results::year_date::YearDate;
use crate::runs::manager::RunManager;
use crate::runs::record::RunStatus;
use eyre::Report;
use itertools::Itertools;
use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use std::collections::BTreeMap;
use std::path::Path;
use treetime_utils::make_report;

/// Request to find a clade of one run in the other time-tree runs.
#[derive(Clone, Debug, Serialize, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CladeRequest {
  /// Run whose tree holds the clade.
  pub run: JobId,
  /// Name of the node at the top of the clade.
  pub node: String,
}

/// A clade of one run found in the other finished time-tree runs.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct CladeInRuns {
  /// Runs whose tree has a node with the same set of samples below it.
  pub matches: Vec<CladeMatch>,
  /// Number of other finished time-tree runs searched.
  pub searched_runs: usize,
  /// Runs whose tree could not be read, with the reason.
  pub unreadable_runs: Vec<UnreadableRun>,
}

/// The node of another run with the same set of samples below it.
#[skip_serializing_none]
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize, JsonSchema)]
pub struct CladeMatch {
  /// Id of the run.
  pub run: JobId,
  /// Title of the run.
  pub title: String,
  /// Name of the node in that run.
  pub node: String,
  /// Date of the node.
  pub date: Option<YearDate>,
  /// Confidence interval of the date.
  pub date_interval: Option<DateInterval>,
}

/// A run whose tree could not be read.
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize, JsonSchema)]
pub struct UnreadableRun {
  /// Id of the run.
  pub run: JobId,
  /// Why the tree could not be read.
  pub message: String,
}

#[derive(Default)]
pub struct TipIds(BTreeMap<String, usize>);

impl TipIds {
  fn id(&mut self, name: &str) -> usize {
    let next = self.0.len();
    *self.0.entry(name.to_owned()).or_insert(next)
  }
}

pub fn clade_keys(tree: &ResultTree, ids: &mut TipIds) -> Vec<Vec<usize>> {
  let mut keys = vec![vec![]; tree.nodes.len()];
  for index in (0..tree.nodes.len()).rev() {
    let node = &tree.nodes[index];
    let mut key = if node.is_tip() {
      vec![ids.id(&node.name)]
    } else {
      node
        .children
        .iter()
        .flat_map(|&child| keys[child].iter().copied())
        .collect()
    };
    key.sort_unstable();
    keys[index] = key;
  }
  keys
}

pub fn first_node_by_key(keys: &[Vec<usize>]) -> BTreeMap<&[usize], usize> {
  let mut first = BTreeMap::new();
  for (index, key) in keys.iter().enumerate() {
    first.entry(key.as_slice()).or_insert(index);
  }
  first
}

pub fn matched_ancestors(first: &ResultTree, second: &ResultTree) -> Vec<(usize, usize)> {
  let mut ids = TipIds::default();
  let first_keys = clade_keys(first, &mut ids);
  let second_keys = clade_keys(second, &mut ids);
  let second_nodes = first_node_by_key(&second_keys);
  first_node_by_key(&first_keys)
    .into_iter()
    .filter(|&(_, index)| !first.nodes[index].is_tip())
    .filter_map(|(key, index)| {
      let other = *second_nodes.get(key)?;
      (!second.nodes[other].is_tip()).then_some((index, other))
    })
    .sorted()
    .collect()
}

pub fn clade_in_trees<'a>(
  tree: &ResultTree,
  node: usize,
  others: impl IntoIterator<Item = &'a ResultTree>,
) -> Vec<Option<usize>> {
  let mut ids = TipIds::default();
  let keys = clade_keys(tree, &mut ids);
  others
    .into_iter()
    .map(|other| {
      let other_keys = clade_keys(other, &mut ids);
      first_node_by_key(&other_keys).get(keys[node].as_slice()).copied()
    })
    .collect()
}

pub fn clade_in_runs(manager: &RunManager, request: &CladeRequest) -> Result<CladeInRuns, Report> {
  let record = manager.get(&request.run)?;
  let tree = read_result_tree(&manager.store().out_dir(&request.run), &record.output_files)?
    .ok_or_else(|| make_report!("run `{}` wrote no Auspice tree", request.run.as_str()))?;
  let node = tree.find(&request.node).ok_or_else(|| {
    make_report!(
      "the tree of run `{}` has no node `{}`",
      request.run.as_str(),
      request.node
    )
  })?;

  let others = manager
    .list()?
    .runs
    .into_iter()
    .filter(|run| run.id != request.run && run.command == AppCommand::Timetree && run.status == RunStatus::Ok)
    .collect::<Vec<_>>();
  let mut unreadable_runs = vec![];
  let mut trees = vec![];
  for run in &others {
    let tree = manager
      .get(&run.id)
      .and_then(|other| read_result_tree(&manager.store().out_dir(&run.id), &other.output_files));
    match tree {
      Ok(Some(tree)) => trees.push((run, tree)),
      Ok(None) => {},
      Err(err) => unreadable_runs.push(UnreadableRun {
        run: run.id.clone(),
        message: format!("{err:#}"),
      }),
    }
  }

  let found = clade_in_trees(&tree, node, trees.iter().map(|(_, tree)| tree));
  let matches = trees
    .iter()
    .zip(found)
    .filter_map(|((run, tree), index)| {
      let other = &tree.nodes[index?];
      Some(CladeMatch {
        run: run.id.clone(),
        title: run.title.clone(),
        node: other.name.clone(),
        date: other.date.clone(),
        date_interval: other.date_interval.clone(),
      })
    })
    .collect();
  Ok(CladeInRuns {
    matches,
    searched_runs: others.len(),
    unreadable_runs,
  })
}

fn read_result_tree(out_dir: &Path, output_files: &[OutputFile]) -> Result<Option<ResultTree>, Report> {
  read_auspice(out_dir, output_files)?
    .map(|auspice| ResultTree::from_auspice(&auspice))
    .transpose()
}
