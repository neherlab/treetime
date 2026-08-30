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
/// partitions (`update_marginal` and `ancestral_reconstruction_marginal` iterate partitions in a
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
  graph: &GraphAncestral,
  index: usize,
  plan: PartitionPlan,
  params: &MarginalPartitionParams,
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

  let sequences = complete_alignment_for_leaves(graph, sequences, &alphabet, params.ignore_missing_alns)?;
  let created = create_marginal_partition(graph, index, alphabet.clone(), &sequences, gtr_model, params.dense)?;
  let partition: Arc<RwLock<dyn MarginalAugurPartition>> = match created.partition {
    MarginalPartition::Sparse(partition) => Arc::new(RwLock::new(partition)),
    MarginalPartition::Dense(partition) => Arc::new(RwLock::new(partition)),
  };

  // Dense partitions attach their leaf sequences here; sparse partitions already carry them from
  // construction (`attach_sequences` is a no-op for sparse), so the call is uniform and safe.
  partition.write_arc().attach_sequences(graph, &sequences)?;

  let single = std::slice::from_ref(&partition);
  update_marginal(graph, single)?;
  ancestral_reconstruction_marginal(
    graph,
    params.include_leaves,
    params.impute_missing_data,
    single,
    params.sample_from_profile,
    rng,
    |_node: &NodeAncestral, _seq: &Seq| Ok(()),
  )?;

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
  pub partition: Arc<RwLock<dyn MarginalAugurPartition>>,
  pub alphabet: Alphabet,
  pub model_name: GtrModelName,
  pub annotation: Option<AugurNodeDataJsonAnnotationEntry>,
  pub reference_override: Option<Seq>,
}

/// A marginal partition that can both take part in the marginal traversal and be read back into augur
/// node data. Both `PartitionMarginalSparse` and `PartitionMarginalDense` satisfy it, so a partition
/// can be erased to `dyn MarginalAugurPartition`.
pub trait MarginalAugurPartition:
  PartitionMarginalOps<NodeAncestral, EdgeAncestral> + AugurNodeDataJsonAncestralPartition
{
}

impl<T> MarginalAugurPartition for T where
  T: PartitionMarginalOps<NodeAncestral, EdgeAncestral> + AugurNodeDataJsonAncestralPartition
{
}
