use crate::ambiguity_outputs::AmbiguityOutputs;
use crate::output_plan::OutputSelection;
use eyre::Report;
use std::collections::{BTreeMap, BTreeSet};
use treetime::alphabet::alphabet::Alphabet;
use treetime::alphabet::state_alphabet::StateAlphabet;
use treetime::ancestral::aa::AaNodeData;
use treetime::partition::marginal::discrete::state_posterior::StatePosterior;
use treetime::partition::marginal::posterior_states::AmbiguityThreshold;
use treetime::seq::mutation::{Mutation, PositionBase};
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::node_names::NodeNames;
use treetime_graph::schedule::{GraphSchedule, TreeSchedule};
use treetime_primitives::Seq;
use treetime_primitives::date::DateConstraint;
use treetime_utils::make_internal_report;
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub struct AnnotatedTreeView<'a> {
  graph: &'a AnnotatedGraph<'a>,
  tree: TreeSchedule<'a>,
}

impl<'a> AnnotatedTreeView<'a> {
  pub fn new(graph: &'a AnnotatedGraph<'a>) -> Result<Self, Report> {
    Ok(Self {
      tree: GraphSchedule::new(graph.graph)?.into_tree()?,
      graph,
    })
  }

  pub fn graph(&self) -> &'a AnnotatedGraph<'a> {
    self.graph
  }

  pub fn tree(&self) -> &TreeSchedule<'a> {
    &self.tree
  }
}

#[derive(Clone, Copy)]
pub struct AnnotatedGraph<'a> {
  pub graph: &'a Graph,
  pub names: &'a NodeNames,
  pub divergence_branch_lengths: &'a BTreeMap<GraphEdgeKey, Option<f64>>,
  pub time_branch_lengths: Option<&'a BTreeMap<GraphEdgeKey, Option<f64>>>,
  pub divergence: Divergence<'a>,
  pub sequences: Option<TreeSequences<'a>>,
  pub dates: Option<TreeDates<'a>>,
  pub traits: Option<TreeTraits<'a>>,
}

impl AnnotatedGraph<'_> {
  pub fn tree_branch_lengths(&self) -> &BTreeMap<GraphEdgeKey, Option<f64>> {
    self.time_branch_lengths.unwrap_or(self.divergence_branch_lengths)
  }
}

#[derive(Clone, Copy)]
pub enum Divergence<'a> {
  CumulativeBranchLength,
  Values(&'a BTreeMap<GraphNodeKey, f64>),
}

#[derive(Clone, Copy)]
pub struct TreeSequences<'a> {
  pub alphabet: &'a Alphabet,
  pub root_sequence: &'a Seq,
  pub edge_mutations: &'a BTreeMap<GraphEdgeKey, Vec<Mutation>>,
  pub amino_acids: Option<TreeAminoAcids<'a>>,
  pub comment_positions: PositionBase,
}

#[derive(Clone, Copy)]
pub struct TreeAminoAcids<'a> {
  pub node_data: &'a AaNodeData,
  pub cdses: &'a BTreeMap<String, AugurNodeDataJsonAnnotationEntry>,
}

#[derive(Clone, Copy)]
pub struct TreeDates<'a> {
  pub num_date: &'a BTreeMap<GraphNodeKey, Option<f64>>,
  pub confidence: Option<&'a BTreeMap<GraphNodeKey, [f64; 2]>>,
  pub excluded: &'a BTreeSet<GraphNodeKey>,
  pub input_dates: Option<&'a BTreeMap<GraphNodeKey, DateConstraint>>,
}

#[derive(Clone, Copy)]
pub struct TreeTraits<'a> {
  pub attribute: &'a str,
  pub states: &'a StateAlphabet,
  pub posteriors: &'a BTreeMap<GraphNodeKey, StatePosterior>,
  pub ambiguity: Option<&'a AmbiguityOutputs>,
}

impl<'a> TreeTraits<'a> {
  pub(crate) fn state_name(&self, key: GraphNodeKey) -> Result<&'a str, Report> {
    let state = self.posteriors[&key].state();
    self
      .states
      .name(state)
      .ok_or_else(|| make_internal_report!("Most likely state {state} of node {key} is not a trait state"))
  }

  pub(crate) fn ambiguity_key(&self, output: OutputSelection) -> Option<String> {
    self
      .ambiguity_threshold(output)
      .map(|_| format!("{}_ambiguous", self.attribute))
  }

  pub(crate) fn ambiguous_states(&self, key: GraphNodeKey, output: OutputSelection) -> Result<Option<String>, Report> {
    let Some(threshold) = self.ambiguity_threshold(output) else {
      return Ok(None);
    };
    let names = threshold
      .plausible_states(self.posteriors[&key].profile())
      .map(|state| {
        self
          .states
          .name(state)
          .ok_or_else(|| make_internal_report!("Plausible state {state} of node {key} is not a trait state"))
      })
      .collect::<Result<Vec<_>, Report>>()?;
    Ok((names.len() > 1).then(|| names.join("|")))
  }

  fn ambiguity_threshold(&self, output: OutputSelection) -> Option<AmbiguityThreshold> {
    self.ambiguity.and_then(|ambiguity| ambiguity.threshold_for(output))
  }
}
