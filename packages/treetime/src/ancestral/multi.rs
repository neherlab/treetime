use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::attach::complete_alignment_for_leaves;
use crate::ancestral::marginal::{ancestral_reconstruction_marginal, update_marginal};
use crate::ancestral::sample::SampleMode;
use crate::gtr::get_gtr::GtrModelName;
use crate::partition::augur::AugurNodeDataJsonAncestralPartition;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::traits::PartitionMarginalOps;
use crate::payload::ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use eyre::Report;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_io::fasta::FastaRecord;
use treetime_primitives::Seq;
use treetime_utils::sync::random::get_random_number_generator;
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
  pub reconstruct_tip_states: bool,
  pub sample_from_profile: SampleMode,
  pub seed: Option<u64>,
  pub ignore_missing_alns: bool,
}

/// Reconstruct several marginal partitions on one shared tree in a single message-passing traversal.
///
/// Multi-partition execution model: every partition is attached to the same graph
/// and the backward+forward marginal passes run once over the whole partition vector, so the tree is
/// walked once rather than once per partition. Per-partition results are then read back by node key
/// (via [`AugurNodeDataJsonAncestralPartition`]), with no cross-graph join by node name.
///
/// Inference is alphabet-agnostic and runs once per partition during construction
/// (`create_marginal_partition` with `--model infer`), matching augur's single `infer_gtr=True`
/// inference; there is no outer GTR-refinement loop here.
pub fn reconstruct_marginal_partitions(
  graph: &GraphAncestral,
  plans: Vec<PartitionPlan>,
  params: &MarginalPartitionParams,
) -> Result<Vec<ReconstructedPartition>, Report> {
  if plans.is_empty() {
    return Ok(Vec::new());
  }

  let mut rng = get_random_number_generator(params.seed);
  let mut partitions: Vec<Arc<RwLock<dyn MarginalAugurPartition>>> = Vec::with_capacity(plans.len());
  let mut metas: Vec<PartitionMeta> = Vec::with_capacity(plans.len());

  for (index, plan) in plans.into_iter().enumerate() {
    let PartitionPlan {
      name,
      alphabet,
      gtr_model,
      sequences,
      annotation,
      reference_override,
    } = plan;

    let sequences = complete_alignment_for_leaves(graph, sequences, &alphabet, params.ignore_missing_alns)?;
    let created = create_marginal_partition(graph, index, alphabet.clone(), &sequences, gtr_model, params.dense)?;
    let partition: Arc<RwLock<dyn MarginalAugurPartition>> = match created.partition {
      MarginalPartition::Sparse(partition) => Arc::new(RwLock::new(partition)),
      MarginalPartition::Dense(partition) => Arc::new(RwLock::new(partition)),
    };

    // Dense partitions attach their leaf sequences here; sparse partitions already carry them from
    // construction (`attach_sequences` is a no-op for sparse), so the call is uniform and safe.
    partition.write_arc().attach_sequences(graph, &sequences)?;

    partitions.push(partition);
    metas.push(PartitionMeta {
      name,
      alphabet,
      model_name: created.model_name,
      annotation,
      reference_override,
    });
  }

  // Single shared backward+forward marginal traversal over all partitions.
  update_marginal(graph, &partitions)?;

  // Resolve the per-node reconstructed states. Reconstruction reads `partitions[0]`, so resolve one
  // partition at a time with a single-element slice; the expensive message passing above was shared.
  for partition in &partitions {
    ancestral_reconstruction_marginal(
      graph,
      params.reconstruct_tip_states,
      std::slice::from_ref(partition),
      params.sample_from_profile,
      &mut rng,
      |_node: &NodeAncestral, _seq: &Seq| Ok(()),
    )?;
  }

  Ok(
    partitions
      .into_iter()
      .zip(metas)
      .map(|(partition, meta)| ReconstructedPartition {
        name: meta.name,
        partition,
        alphabet: meta.alphabet,
        model_name: meta.model_name,
        annotation: meta.annotation,
        reference_override: meta.reference_override,
      })
      .collect(),
  )
}

/// A reconstructed partition and the metadata needed to serialize it into augur node data.
pub struct ReconstructedPartition {
  pub name: String,
  pub partition: Arc<RwLock<dyn MarginalAugurPartition>>,
  pub alphabet: Alphabet,
  pub model_name: GtrModelName,
  pub annotation: Option<AugurNodeDataJsonAnnotationEntry>,
  pub reference_override: Option<Seq>,
}

/// A marginal partition that can both take part in the marginal traversal and be read back into augur
/// node data. Both `PartitionMarginalSparse` and `PartitionMarginalDense` satisfy it, so a
/// heterogeneous partition vector can be erased to `dyn MarginalAugurPartition`.
pub trait MarginalAugurPartition:
  PartitionMarginalOps<NodeAncestral, EdgeAncestral> + AugurNodeDataJsonAncestralPartition
{
}

impl<T> MarginalAugurPartition for T where
  T: PartitionMarginalOps<NodeAncestral, EdgeAncestral> + AugurNodeDataJsonAncestralPartition
{
}

struct PartitionMeta {
  name: String,
  alphabet: Alphabet,
  model_name: GtrModelName,
  annotation: Option<AugurNodeDataJsonAnnotationEntry>,
  reference_override: Option<Seq>,
}
