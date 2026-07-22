use crate::gtr::gtr::GTR;
use crate::gtr::infer_gtr::common::MutationCounts;
use crate::make_error;
use crate::partition::dense::{DenseEdgePartition, DenseNodePartition, DenseSeqDistribution, DenseSeqInfo};
use crate::partition::discrete_states::DiscreteStates;
use crate::partition::marginal_core::{
  IndexedMarginalPartition, MarginalData, MarginalPartition, marginal_process_backward_indexed,
  marginal_process_forward_indexed,
};
use crate::partition::traits::{HasGtr, HasLogLh, PartitionMarginalPasses, TransitionCounting};
use eyre::Report;
use indexmap::IndexSet;
use itertools::Itertools;
use log::warn;
use maplit::btreemap;
use ndarray::{Array1, Array2};
use serde::Serialize;
use std::collections::BTreeMap;
use treetime_graph::edge::EdgeOptimizeOps;
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, GraphNodeKey, Named};
use treetime_io::nwk::NodeCommentProvider;
use treetime_utils::array::ndarray::argmax_first;

#[derive(Clone, Debug, Serialize)]
pub struct PartitionMarginalDiscrete {
  pub data: MarginalData,
  pub states: DiscreteStates,
}

impl PartitionMarginalDiscrete {
  pub fn new(gtr: GTR, states: DiscreteStates, min_branch_length: f64, filter_uninformative_root: bool) -> Self {
    Self {
      data: MarginalData {
        gtr,
        nodes: btreemap! {},
        edges: btreemap! {},
        min_branch_length,
        filter_uninformative_root,
      },
      states,
    }
  }

  pub fn n_states(&self) -> usize {
    self.states.len()
  }

  pub fn attach_traits<N, E>(
    &mut self,
    graph: &Graph<N, E, ()>,
    traits: &BTreeMap<String, String>,
  ) -> Result<(), Report>
  where
    N: GraphNode + Named,
    E: EdgeOptimizeOps,
  {
    let n_states = self.n_states();
    validate_trait_names(graph, traits)?;

    for leaf in graph.get_leaves() {
      let leaf_guard = leaf.read_arc();
      let leaf_key = leaf_guard.key();
      let leaf_payload = leaf_guard.payload().read_arc();
      let leaf_name = leaf_payload.name().map(|n| n.as_ref().to_owned()).unwrap_or_default();

      let profile = if let Some(trait_value) = traits.get(&leaf_name) {
        if let Some(index) = self.states.get_index(trait_value) {
          one_hot_profile(index, n_states)
        } else {
          uniform_profile(n_states)
        }
      } else {
        uniform_profile(n_states)
      };

      self.data.nodes.insert(
        leaf_key,
        DenseNodePartition {
          seq: DenseSeqInfo::default(),
          profile: DenseSeqDistribution::new(profile, 0.0),
        },
      );
    }

    for edge in graph.get_edges() {
      let edge_key = edge.read_arc().key();
      self.data.edges.insert(edge_key, DenseEdgePartition::default());
    }

    Ok(())
  }

  pub fn get_reconstructed_trait(&self, node_key: GraphNodeKey) -> Option<String> {
    let node = self.data.nodes.get(&node_key)?;
    let row = node.profile.dis.row(0);
    let argmax = argmax_first(&row)?;
    Some(self.states.get_name(argmax).to_owned())
  }

  pub fn get_confidence(&self, node_key: GraphNodeKey) -> Option<Array1<f64>> {
    let node = self.data.nodes.get(&node_key)?;
    Some(node.profile.dis.row(0).to_owned())
  }
}

impl HasGtr for PartitionMarginalDiscrete {
  fn gtr(&self) -> &GTR {
    &self.data.gtr
  }

  fn gtr_mut(&mut self) -> &mut GTR {
    &mut self.data.gtr
  }

  fn sequence_length(&self) -> usize {
    1
  }
}

impl HasLogLh for PartitionMarginalDiscrete {
  fn get_log_lh(&self, node_key: GraphNodeKey) -> f64 {
    self.data.nodes.get(&node_key).map_or(0.0, |node| node.profile.log_lh)
  }

  fn reset_node_log_likelihoods(&mut self) {
    for node_data in self.data.nodes.values_mut() {
      node_data.profile.log_lh = 0.0;
    }
  }
}

impl<N, E> TransitionCounting<N, E> for PartitionMarginalDiscrete
where
  N: GraphNode,
  E: EdgeOptimizeOps,
{
  fn count_transitions(&self, graph: &Graph<N, E, ()>) -> Result<MutationCounts, Report> {
    self.data.count_transitions(graph)
  }
}

impl<N, E> MarginalPartition<N, E> for PartitionMarginalDiscrete
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn marginal_data(&self) -> &MarginalData {
    &self.data
  }

  fn marginal_data_mut(&mut self) -> &mut MarginalData {
    &mut self.data
  }

  fn indexed_storage_mut(
    &mut self,
  ) -> (
    &mut BTreeMap<GraphNodeKey, DenseNodePartition>,
    &mut BTreeMap<treetime_graph::edge::GraphEdgeKey, DenseEdgePartition>,
  ) {
    (&mut self.data.nodes, &mut self.data.edges)
  }

  fn leaf_profile(&self, node_key: GraphNodeKey) -> Result<DenseSeqDistribution, Report> {
    let node = &self.data.nodes[&node_key];
    Ok(node.profile.clone())
  }
}

impl<N, E> IndexedMarginalPartition<N, E> for PartitionMarginalDiscrete
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn indexed_missing_node(&self, _key: GraphNodeKey) -> Result<DenseNodePartition, Report> {
    Ok(DenseNodePartition {
      seq: DenseSeqInfo::default(),
      profile: DenseSeqDistribution::default(),
    })
  }

  fn indexed_leaf_profile(&self, node: &DenseNodePartition) -> Result<DenseSeqDistribution, Report> {
    Ok(node.profile.clone())
  }

  fn indexed_forward_post(
    &self,
    _is_root: bool,
    _is_leaf: bool,
    _parent: Option<&DenseNodePartition>,
    _node: &mut DenseNodePartition,
    _parent_edge: Option<&mut DenseEdgePartition>,
  ) -> Result<(), Report> {
    Ok(())
  }
}

impl<N, E> PartitionMarginalPasses<N, E> for PartitionMarginalDiscrete
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  fn process_backward_pass(&mut self, graph: &Graph<N, E, ()>) -> Result<(), Report> {
    marginal_process_backward_indexed(self, graph)
  }

  fn process_forward_pass(&mut self, graph: &Graph<N, E, ()>) -> Result<(), Report> {
    marginal_process_forward_indexed(self, graph)
  }

  fn get_sequence_length(&self) -> usize {
    1
  }
}

pub struct DiscreteCommentProvider<'a> {
  partition: &'a PartitionMarginalDiscrete,
  attribute: &'a str,
}

impl<'a> DiscreteCommentProvider<'a> {
  pub fn new(partition: &'a PartitionMarginalDiscrete, attribute: &'a str) -> Self {
    Self { partition, attribute }
  }
}

impl NodeCommentProvider for DiscreteCommentProvider<'_> {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    Ok(
      self
        .partition
        .get_reconstructed_trait(key)
        .map_or_else(BTreeMap::new, |trait_value| {
          btreemap! {
            self.attribute.to_owned() => trait_value,
          }
        }),
    )
  }
}

fn one_hot_profile(index: usize, n_states: usize) -> Array2<f64> {
  let mut profile = Array2::zeros((1, n_states));
  profile[[0, index]] = 1.0;
  profile
}

fn uniform_profile(n_states: usize) -> Array2<f64> {
  Array2::from_elem((1, n_states), 1.0 / n_states as f64)
}

fn validate_trait_names<N, E>(graph: &Graph<N, E, ()>, traits: &BTreeMap<String, String>) -> Result<(), Report>
where
  N: GraphNode + Named,
  E: EdgeOptimizeOps,
{
  let leaf_names: IndexSet<String> = graph
    .get_leaves()
    .iter()
    .map(|leaf| {
      let leaf = leaf.read_arc();
      let payload = leaf.payload().read_arc();
      payload.name().map(|name| name.as_ref().to_owned()).unwrap_or_default()
    })
    .collect();
  let trait_names: IndexSet<String> = traits.keys().cloned().collect();

  let missing_in_metadata: IndexSet<String> = leaf_names.difference(&trait_names).cloned().collect();
  if !missing_in_metadata.is_empty() {
    return make_error!(
      "Mugration: tree leaves missing from metadata: {}",
      missing_in_metadata.iter().join(", ")
    );
  }

  // Metadata naming samples absent from the tree is the common case: metadata files are shared
  // across analyses and routinely list more samples than a pruned or subsampled tree. Warn for
  // visibility (matching the dates subsystem) instead of rejecting the run.
  let missing_in_tree: IndexSet<String> = trait_names.difference(&leaf_names).cloned().collect();
  if !missing_in_tree.is_empty() {
    let sample = missing_in_tree.iter().take(10).join(", ");
    let suffix = if missing_in_tree.len() > 10 { "..." } else { "" };
    warn!(
      "Mugration: {} metadata names not present in tree: {sample}{suffix}",
      missing_in_tree.len()
    );
  }

  Ok(())
}
