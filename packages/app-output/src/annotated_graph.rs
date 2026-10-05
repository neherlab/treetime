use eyre::Report;
use ndarray::Array1;
use std::collections::{BTreeMap, BTreeSet};
use treetime::ancestral::aa::AaNodeData;
use treetime::partition::storage::discrete::DiscreteStates;
use treetime::seq::mutation::Mutation;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::tree_view::TreeView;
use treetime_primitives::Seq;
use treetime_primitives::date::DatesMap;
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub struct AnnotatedGraph<'a> {
  pub graph: &'a Graph,
  pub names: &'a BTreeMap<GraphNodeKey, Option<String>>,
  pub divergence_branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  pub time_branch_lengths: Option<&'a BTreeMap<GraphEdgeKey, Option<f64>>>,
  pub divergence: Divergence<'a>,
  pub branch_support: Option<&'a BTreeMap<GraphNodeKey, Option<f64>>>,
  pub sequences: Option<TreeSequences<'a>>,
  pub dates: Option<TreeDates<'a>>,
  pub traits: Option<TreeTraits<'a>>,
}

impl AnnotatedGraph<'_> {
  pub fn tree_branch_lengths(&self) -> &BTreeMap<GraphEdgeKey, Option<f64>> {
    self.time_branch_lengths.unwrap_or(self.divergence_branch_lengths)
  }
}

pub struct AnnotatedTreeView<'a> {
  graph: &'a AnnotatedGraph<'a>,
  tree: TreeView<'a>,
}

impl<'a> AnnotatedTreeView<'a> {
  pub fn new(graph: &'a AnnotatedGraph<'a>) -> Result<Self, Report> {
    Ok(Self {
      tree: TreeView::new(graph.graph)?,
      graph,
    })
  }

  pub fn graph(&self) -> &'a AnnotatedGraph<'a> {
    self.graph
  }

  pub fn tree(&self) -> &TreeView<'a> {
    &self.tree
  }
}

pub enum Divergence<'a> {
  CumulativeBranchLength,
  Values(&'a BTreeMap<GraphNodeKey, f64>),
}

pub struct TreeSequences<'a> {
  pub root_sequence: &'a Seq,
  pub edge_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  pub mutation_counts: Option<&'a BTreeMap<GraphEdgeKey, usize>>,
  pub amino_acids: Option<TreeAminoAcids<'a>>,
}

pub struct TreeAminoAcids<'a> {
  pub node_data: &'a AaNodeData,
  pub cdses: &'a BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
}

pub struct TreeDates<'a> {
  pub num_date: &'a BTreeMap<GraphNodeKey, Option<f64>>,
  pub confidence: Option<&'a BTreeMap<GraphNodeKey, [f64; 2]>>,
  pub excluded: &'a BTreeSet<GraphNodeKey>,
  pub input_dates: Option<&'a DatesMap>,
}

pub struct TreeTraits<'a> {
  pub attribute: &'a str,
  pub states: &'a DiscreteStates,
  pub values: &'a BTreeMap<GraphNodeKey, Option<String>>,
  pub profiles: &'a BTreeMap<GraphNodeKey, Option<Array1<f64>>>,
}
