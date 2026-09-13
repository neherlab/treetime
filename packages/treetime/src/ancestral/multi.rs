use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::attach::complete_alignment_for_leaves;
use crate::ancestral::marginal::{ancestral_reconstruction, profile_branch_lengths};
use crate::ancestral::pipeline::{DenseReconstruction, SparseReconstruction};
use crate::ancestral::sample::SampleMode;
use crate::gtr::get_gtr::GtrModelName;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::io::augur::AugurNodeDataJsonAncestralPartition;
use crate::partition::marginal::shared::update::MarginalUpdate;
use eyre::Report;
use std::collections::BTreeMap;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaRecord;
use treetime_primitives::Seq;
use util_augur_node_data_json::AugurNodeDataJsonAnnotationEntry;

/// A partition to reconstruct on the shared tree.
///
/// `name` is the key under which the partition appears in the augur node-data JSON (`nuc` for the
/// nucleotide partition, the CDS name for amino-acid partitions). `sequences` are the per-leaf
/// sequences for this partition only; different partitions can have different lengths and alphabets.
pub struct PartitionPlan {
  pub name: String,
  pub alphabet: Alphabet,
  pub gtr_model: GtrModelName,
  pub sequences: Vec<FastaRecord>,
  pub annotation: Option<AugurNodeDataJsonAnnotationEntry>,
  pub reference_override: Option<Seq>,
}

/// Reconstruction parameters shared by every partition in a multi-partition run.
pub struct MarginalPartitionParams {
  pub dense: Option<bool>,
  pub include_leaves: bool,
  pub impute_missing_data: bool,
  pub sample_from_profile: SampleMode,
  pub seed: Option<u64>,
  pub ignore_missing_alns: bool,
}

/// Reconstruct one marginal partition on the shared tree and return its per-node results.
///
/// Each partition is independent: its GTR inference, backward and forward marginal passes, and node
/// reconstruction touch only this partition's own message state, with no data shared across
/// partitions (`marginal_update` and `ancestral_reconstruction_marginal` iterate partitions in a
/// plain loop, so a single-element slice does the same work as a slice of many). Callers that
/// reconstruct several partitions (per-CDS amino-acid alignments) therefore call this once per
/// partition and consume each result before building the next, so the resident marginal state stays
/// bounded to a single partition instead of scaling with the partition count.
///
/// `index` distinguishes partitions during construction. `rng` is passed in by the caller so that
/// sampled reconstruction (`--sample-from-profile=root|all`) draws in a fixed partition order.
///
/// Inference is alphabet-agnostic and runs once during construction (`create_marginal_partition`
/// with `--model infer`), matching augur's single `infer_gtr=True` inference; there is no outer
/// GTR-refinement loop here.
pub fn reconstruct_marginal_partition(
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
  let created = create_marginal_partition(
    graph,
    index,
    alphabet.clone(),
    &sequences,
    gtr_model,
    params.dense,
    branch_lengths,
    names,
  )?;
  let profile_lengths = profile_branch_lengths(branch_lengths);

  // Each partition runs its own marginal passes and node reconstruction over its own role-typed result
  // maps, then is boxed as the augur read view. The two representations reconstruct differently (sparse
  // tip imputation reads the forward down-message), so each branch owns its reconstruction closure.
  let partition: Box<dyn AugurNodeDataJsonAncestralPartition> = match created.partition {
    MarginalPartition::Sparse(partition, node_states) => {
      let MarginalUpdate {
        mut node_states,
        backward,
        forward,
        estimates,
        log_lh: _,
      } = partition.marginal_update(graph, &profile_lengths, node_states)?;
      ancestral_reconstruction(
        graph,
        |node| {
          partition.reconstruct_node_sequence(
            &mut node_states,
            &forward,
            node,
            params.include_leaves,
            params.impute_missing_data,
            params.sample_from_profile,
            rng,
          )
        },
        |_key: GraphNodeKey, _seq: &Seq| Ok(()),
      )?;
      Box::new(SparseReconstruction {
        partition,
        node_states,
        backward,
        forward,
        estimates,
      })
    },
    MarginalPartition::Dense(partition) => {
      let node_states = partition.attach_sequences(graph, &sequences, names)?;
      let MarginalUpdate {
        mut node_states,
        backward,
        forward,
        estimates,
        log_lh: _,
      } = partition.marginal_update(graph, &profile_lengths, node_states)?;
      ancestral_reconstruction(
        graph,
        |node| {
          partition.reconstruct_node_sequence(
            &mut node_states,
            node,
            params.include_leaves,
            params.impute_missing_data,
            params.sample_from_profile,
            rng,
          )
        },
        |_key: GraphNodeKey, _seq: &Seq| Ok(()),
      )?;
      Box::new(DenseReconstruction {
        partition,
        node_states,
        backward,
        forward,
        estimates,
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

/// A reconstructed partition and the metadata needed to serialize it into augur node data.
pub struct ReconstructedPartition {
  pub name: String,
  pub partition: Box<dyn AugurNodeDataJsonAncestralPartition>,
  pub alphabet: Alphabet,
  pub model_name: GtrModelName,
  pub annotation: Option<AugurNodeDataJsonAnnotationEntry>,
  pub reference_override: Option<Seq>,
}
