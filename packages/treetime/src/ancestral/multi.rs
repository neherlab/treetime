use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::attach::complete_alignment_for_leaves;
use crate::ancestral::marginal::{ancestral_reconstruction, branch_lengths_or_zero};
use crate::ancestral::pipeline::{AncestralPartition, DenseReconstruction, SparseReconstruction};
use crate::ancestral::sample::SampleMode;
use crate::ancestral::tip_states::TipStates;
use crate::gtr::get_gtr::GtrModelName;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::marginal::shared::update::{MarginalPasses, MarginalUpdate};
use crate::seq::alignment::node_seq_inputs;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_primitives::{AlignmentRecord, Seq};
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

pub(crate) fn reconstruct_marginal_partition(
  graph: &Graph,
  index: usize,
  plan: PartitionPlan,
  params: &MarginalPartitionParams,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  rng: &mut dyn rand::RngCore,
) -> Result<ReconstructedPartition, Report> {
  let PartitionPlan {
    name,
    alphabet,
    gtr_model,
    sequences,
    annotation,
    reference_override,
  } = plan;

  let sequences = complete_alignment_for_leaves(graph, sequences, &alphabet, params.ignore_missing_alns, names)?;
  let node_inputs = node_seq_inputs(graph, names, sequences);
  let profile_lengths = branch_lengths_or_zero(branch_lengths);
  let created = create_marginal_partition(
    graph,
    index,
    alphabet.clone(),
    &node_inputs,
    gtr_model,
    params.dense,
    &profile_lengths,
  )?;
  let gtr = created.gtr;

  let partition: AncestralPartition = match created.partition {
    MarginalPartition::Sparse(partition, node_states) => {
      let MarginalUpdate {
        mut node_states, edges, ..
      } = partition.marginal_update(&gtr, graph, &profile_lengths, node_states)?;
      ancestral_reconstruction(graph, |node| {
        partition.advance_node_state(
          &mut node_states,
          &edges.forward,
          node,
          TipStates {
            include_leaves: params.include_leaves,
            impute: params.impute_missing_data,
          },
          params.sample_from_profile,
          rng,
        )
      })?;
      AncestralPartition::Sparse(SparseReconstruction {
        partition,
        gtr,
        node_states,
        edges,
      })
    },
    MarginalPartition::Dense(partition) => {
      let node_states = partition.attach_sequences(graph, &node_inputs)?;
      let MarginalUpdate {
        mut node_states, edges, ..
      } = partition.marginal_update(&gtr, graph, &profile_lengths, node_states)?;
      ancestral_reconstruction(graph, |node| {
        partition
          .reconstruct_node_sequence(
            &mut node_states,
            node,
            TipStates {
              include_leaves: params.include_leaves,
              impute: params.impute_missing_data,
            },
            params.sample_from_profile,
            rng,
          )
          .map(|_| ())
      })?;
      AncestralPartition::Dense(DenseReconstruction {
        partition,
        gtr,
        node_states,
        edges,
      })
    },
  };

  Ok(ReconstructedPartition {
    name,
    partition,
    alphabet,
    model_name: created.model_name,
    annotation,
    reference_override,
  })
}

pub struct PartitionPlan {
  pub name: String,
  pub alphabet: Alphabet,
  pub gtr_model: GtrModelName,
  pub sequences: Vec<AlignmentRecord>,
  pub annotation: Option<AugurNodeDataJsonAnnotationEntry>,
  pub reference_override: Option<Seq>,
}

pub struct MarginalPartitionParams {
  pub dense: Option<bool>,
  pub include_leaves: bool,
  pub impute_missing_data: bool,
  pub sample_from_profile: SampleMode,
  pub seed: Option<u64>,
  pub ignore_missing_alns: bool,
}

pub struct ReconstructedPartition {
  pub name: String,
  pub partition: AncestralPartition,
  pub alphabet: Alphabet,
  pub model_name: GtrModelName,
  pub annotation: Option<AugurNodeDataJsonAnnotationEntry>,
  pub reference_override: Option<Seq>,
}
