use crate::edge::GraphEdgeKey;
use crate::graph::Graph;
use crate::key_map::{EdgeMap, NodeMap};
use crate::node::GraphNodeKey;
use crate::node_names::{InputNodeNames, NodeNames};
use crate::schedule::{GraphSchedule, find_cycle};
use crate::scope::{EdgeScopes, SiteScope};
use eyre::{Report, WrapErr};
use indexmap::IndexMap;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::btree_map::Entry;
use treetime_utils::{make_error, make_internal_report};

pub const SHARED_BRANCH_LENGTH_TOLERANCE: f64 = 1e-6;

pub struct NamedTree<'a> {
  pub label: &'a str,
  pub graph: &'a Graph,
  pub names: &'a NodeNames,
  pub unnamed_nodes: &'a BTreeSet<GraphNodeKey>,
  pub branch_lengths: &'a EdgeMap<Option<f64>>,
}

#[derive(Debug)]
pub struct MergedGraph {
  pub graph: Graph,
  pub names: NodeNames,
  pub branch_lengths: EdgeMap<Option<f64>>,
  pub scopes: EdgeScopes,
}

pub fn merge_trees_by_name(trees: &[NamedTree<'_>]) -> Result<MergedGraph, Report> {
  let mut graph = Graph::new();
  let mut by_name: IndexMap<&str, GraphNodeKey> = IndexMap::new();
  let mut edges: BTreeMap<(GraphNodeKey, GraphNodeKey), GraphEdgeKey> = BTreeMap::new();
  let mut branch_lengths = EdgeMap::new();
  let mut scopes = EdgeScopes::new(trees.len());

  for (segment, tree) in trees.iter().enumerate() {
    let schedule = GraphSchedule::new(tree.graph).wrap_err_with(|| format!("When reading tree '{}'", tree.label))?;
    schedule
      .require_tree()
      .wrap_err_with(|| format!("When reading tree '{}'", tree.label))?;
    require_unique_names(tree)?;

    for &node in schedule.preorder() {
      let merged = *by_name.entry(&tree.names[node]).or_insert_with(|| graph.add_node());
      let Some(&(parent, edge)) = schedule.parents(node).first() else {
        continue;
      };
      let merged_parent = *by_name
        .get(&tree.names[parent])
        .ok_or_else(|| make_internal_report!("Parent {parent} comes after its child {node} in the preorder"))?;
      let length = tree.branch_lengths.get(edge).copied().flatten();
      let merged_edge = match edges.entry((merged_parent, merged)) {
        Entry::Vacant(entry) => {
          let key = graph.add_edge(merged_parent, merged)?;
          branch_lengths.insert(key, length);
          *entry.insert(key)
        },
        Entry::Occupied(entry) => {
          let key = *entry.get();
          let shared = branch_lengths
            .get_mut(key)
            .ok_or_else(|| make_internal_report!("Merged edge {key} has no branch length entry"))?;
          *shared = shared_branch_length(*shared, length).wrap_err_with(|| {
            format!(
              "When merging the edge '{}' -> '{}' of tree '{}' into the same edge of an earlier tree",
              &tree.names[parent], &tree.names[node], tree.label
            )
          })?;
          key
        },
      };
      scopes.set(merged_edge, segment, SiteScope::Full)?;
    }
  }

  graph.build()?;
  let input_names: NodeMap<Option<String>> = by_name
    .iter()
    .map(|(name, &key)| (key, Some((*name).to_owned())))
    .collect();
  let InputNodeNames { names, .. } = NodeNames::from_input(input_names, &graph)?;

  if let Some(cycle) = find_cycle(&graph)? {
    let path = cycle
      .iter()
      .rev()
      .chain(cycle.last())
      .map(|&key| &names[key])
      .join("' -> '");
    return make_error!(
      "The trees order their nodes in conflicting ways: merged by name they form the cycle '{path}', where each node is a parent of the next"
    );
  }
  let schedule = GraphSchedule::new(&graph)?;
  scopes.validate(&schedule, |key| names[key].to_owned())?;

  Ok(MergedGraph {
    graph,
    names,
    branch_lengths,
    scopes,
  })
}

fn require_unique_names(tree: &NamedTree<'_>) -> Result<(), Report> {
  if !tree.unnamed_nodes.is_empty() {
    return make_error!(
      "The tree '{}' has {} nodes without a name. Trees merge by node name, so every node needs one",
      tree.label,
      tree.unnamed_nodes.len()
    );
  }
  let repeated = tree
    .graph
    .node_keys()
    .filter(|&key| tree.names.pairing_name(key) != &tree.names[key])
    .map(|key| tree.names.pairing_name(key))
    .unique()
    .collect_vec();
  if !repeated.is_empty() {
    return make_error!(
      "The tree '{}' gives the same name to more than one node: '{}'. Trees merge by node name, so names must be unique within a tree",
      tree.label,
      repeated.join("', '")
    );
  }
  Ok(())
}

fn shared_branch_length(earlier: Option<f64>, later: Option<f64>) -> Result<Option<f64>, Report> {
  match (earlier, later) {
    (Some(a), Some(b)) => {
      let tolerance = SHARED_BRANCH_LENGTH_TOLERANCE * a.abs().max(b.abs()).max(1.0);
      if (a - b).abs() > tolerance {
        return make_error!("The two trees give the edge different branch lengths, {a} and {b}");
      }
      Ok(Some(a))
    },
    (earlier, later) => Ok(earlier.or(later)),
  }
}
