use crate::edge::GraphEdgeKey;
use crate::graph::Graph;
use crate::key_map::{EdgeMap, NodeMap};
use crate::node::GraphNodeKey;
use crate::node_names::{InputNodeNames, NodeNames};
use crate::pass_through::{PassThroughReport, SegmentTree, insert_pass_through_nodes};
use crate::schedule::{GraphSchedule, find_cycle};
use crate::scope::{EdgeScopes, SiteScope};
use eyre::{Report, WrapErr};
use indexmap::IndexMap;
use itertools::Itertools;
use std::collections::BTreeMap;
use std::collections::BTreeSet;
use std::collections::btree_map::Entry;
use treetime_utils::{make_error, make_internal_report};

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
  pub pass_through: PassThroughReport,
}

pub fn merge_trees_by_name(trees: &[NamedTree<'_>]) -> Result<MergedGraph, Report> {
  let mut segments = trees.iter().map(segment_tree).collect::<Result<Vec<_>, Report>>()?;
  let pass_through = insert_pass_through_nodes(&mut segments);

  let mut graph = Graph::new();
  let mut by_name: IndexMap<&str, GraphNodeKey> = IndexMap::new();
  let mut edges: BTreeMap<(GraphNodeKey, GraphNodeKey), GraphEdgeKey> = BTreeMap::new();
  let mut length_sums: EdgeMap<LengthSum> = EdgeMap::new();
  let mut scopes = EdgeScopes::new(trees.len());

  for (segment, tree) in segments.iter().enumerate() {
    for name in tree.preorder() {
      let merged = *by_name.entry(name).or_insert_with(|| graph.add_node());
      let Some(&(parent, length)) = tree.parent.get(name) else {
        continue;
      };
      let merged_parent = *by_name
        .get(parent)
        .ok_or_else(|| make_internal_report!("Parent '{parent}' comes after its child '{name}' in the preorder"))?;
      let merged_edge = match edges.entry((merged_parent, merged)) {
        Entry::Vacant(entry) => *entry.insert(graph.add_edge(merged_parent, merged)?),
        Entry::Occupied(entry) => *entry.get(),
      };
      if let Some(length) = length {
        let sum = length_sums.entry(merged_edge).or_default();
        sum.total += length;
        sum.count += 1.0;
      }
      scopes.set(merged_edge, segment, SiteScope::Full)?;
    }
  }

  graph.build()?;
  let branch_lengths = graph
    .get_edges()
    .map(|edge| {
      let mean = length_sums.get(edge.key()).map(|sum| sum.total / sum.count);
      (edge.key(), mean)
    })
    .collect();
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
    pass_through,
  })
}

fn segment_tree<'a>(tree: &NamedTree<'a>) -> Result<SegmentTree<'a>, Report> {
  let schedule = GraphSchedule::new(tree.graph).wrap_err_with(|| format!("When reading tree '{}'", tree.label))?;
  schedule
    .require_tree()
    .wrap_err_with(|| format!("When reading tree '{}'", tree.label))?;
  require_unique_names(tree)?;
  let names = tree.names;
  let root = schedule
    .roots()
    .first()
    .map(|&root| &names[root])
    .ok_or_else(|| make_internal_report!("Tree '{}' has no root", tree.label))?;
  let mut children = BTreeMap::new();
  let mut parent = BTreeMap::new();
  for &node in schedule.preorder() {
    let kids = schedule
      .children(node)
      .iter()
      .map(|&(child, _)| &names[child])
      .collect_vec();
    if !kids.is_empty() {
      children.insert(&names[node], kids);
    }
    if let Some(&(up, edge)) = schedule.parents(node).first() {
      parent.insert(
        &names[node],
        (&names[up], tree.branch_lengths.get(edge).copied().flatten()),
      );
    }
  }
  Ok(SegmentTree { root, children, parent })
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

#[derive(Default)]
struct LengthSum {
  total: f64,
  count: f64,
}
