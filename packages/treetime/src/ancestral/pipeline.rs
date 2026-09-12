use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::attach::complete_alignment_for_leaves;
use crate::ancestral::fitch::{ancestral_reconstruction_fitch, create_fitch_partition};
use crate::ancestral::marginal::{ancestral_reconstruction, profile_branch_lengths};
use crate::ancestral::mask::create_mask;
use crate::ancestral::params::MethodAncestral;
use crate::ancestral::sample::SampleMode;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::gtr::refinement::refine_gtr_iterative;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::marginal::dense::partition::PartitionMarginalDense;
use crate::partition::marginal::sparse::partition::PartitionMarginalSparse;
use crate::partition::storage::dense::{DenseEdgeBackward, DenseEdgeEstimate, DenseEdgeForward, DenseNodeState};
use crate::partition::storage::sparse::{SparseEdgeBackward, SparseEdgeForward, SparseNodeState};
use crate::partition::traits::HasGtr;
use crate::progress::ProgressSink;
use crate::seq::alignment::get_common_length;
use crate::seq::mutation::Sub;
use eyre::Report;
use serde::Serialize;
use std::collections::BTreeMap;
use strum::VariantNames;
use treetime_graph::edge::GraphEdgeKey;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaRecord;
use treetime_primitives::Seq;
use treetime_utils::make_error;
use treetime_utils::sync::random::get_random_number_generator;

pub struct AncestralParams {
  pub method: MethodAncestral,
  pub model: GtrModelName,
  pub dense: Option<bool>,
  /// Emit reconstructed leaf sequences in addition to internal nodes.
  pub include_leaves: bool,
  /// Resolve ambiguous and unknown tip states (`N` and IUPAC codes) to the most likely inferred
  /// state. Only defined for marginal reconstruction; a no-op for Fitch parsimony.
  pub impute_missing_data: bool,
  pub gtr_iterations: usize,
  pub site_specific_gtr: bool,
  pub seed: Option<u64>,
  pub sample_from_profile: SampleMode,
  pub ignore_missing_alns: bool,
}

pub struct AncestralInput {
  pub graph: Graph,
  pub alphabet: Alphabet,
  pub sequences: Vec<FastaRecord>,
}

/// A completed sparse reconstruction: the durable partition inputs together with the node states and
/// per-edge messages/estimates the passes returned, kept as distinct owned maps. The output writers
/// build a short-lived read view over these.
#[derive(Clone, Serialize)]
pub struct SparseReconstruction {
  pub partition: PartitionMarginalSparse,
  pub node_states: BTreeMap<GraphNodeKey, SparseNodeState>,
  pub backward: BTreeMap<GraphEdgeKey, SparseEdgeBackward>,
  pub forward: BTreeMap<GraphEdgeKey, SparseEdgeForward>,
  pub estimates: BTreeMap<GraphEdgeKey, Vec<Sub>>,
}

/// A completed dense reconstruction: the durable partition inputs together with the node states and
/// per-edge messages/estimates the passes returned, kept as distinct owned maps.
#[derive(Clone, Serialize)]
pub struct DenseReconstruction {
  pub partition: PartitionMarginalDense,
  pub node_states: BTreeMap<GraphNodeKey, DenseNodeState>,
  pub backward: BTreeMap<GraphEdgeKey, DenseEdgeBackward>,
  pub forward: BTreeMap<GraphEdgeKey, DenseEdgeForward>,
  pub estimates: BTreeMap<GraphEdgeKey, DenseEdgeEstimate>,
}

#[derive(Clone, Serialize)]
#[serde(rename_all = "kebab-case")]
pub enum AncestralPartition {
  Fitch(PartitionFitch),
  Sparse(SparseReconstruction),
  Dense(DenseReconstruction),
}

#[derive(Debug, Serialize)]
pub struct AncestralOutput {
  #[serde(skip)]
  pub graph: Graph,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub mask: Vec<bool>,
  /// Reconstructed sequences keyed by node id, captured from the serial reconstruction walk.
  #[serde(skip)]
  pub node_sequences: BTreeMap<GraphNodeKey, Seq>,
}

pub struct AncestralOutputFull {
  pub output: AncestralOutput,
  pub partition: Option<AncestralPartition>,
}

pub fn run<F>(
  params: &AncestralParams,
  input: AncestralInput,
  names: &BTreeMap<GraphNodeKey, Option<String>>,
  branch_lengths: &BTreeMap<GraphEdgeKey, Option<f64>>,
  mut on_sequence: F,
  progress: &dyn ProgressSink,
) -> Result<AncestralOutputFull, Report>
where
  F: FnMut(GraphNodeKey, &Seq) -> Result<(), Report>,
{
  let profile_lengths = profile_branch_lengths(branch_lengths);
  if params.site_specific_gtr {
    return make_error!(
      "--site-specific-gtr is not yet integrated into the ancestral reconstruction pipeline. \
       The mathematical core (GTRSiteSpecific) is implemented but partition system wiring is pending."
    );
  }

  if params.sample_from_profile != SampleMode::Argmax && params.method != MethodAncestral::Marginal {
    return make_error!(
      "--sample-from-profile={:?} requires --method-anc=marginal. Posterior sampling is only defined \
       for marginal reconstruction; {:?} has no posterior profile to sample.",
      params.sample_from_profile,
      params.method
    );
  }

  let AncestralInput {
    graph,
    alphabet,
    sequences,
  } = input;

  // Tips absent from the alignment become fully-ambiguous sequences here, once, for every
  // partition backend (fitch, sparse, dense) and alphabet (nucleotide, amino acid). After this the
  // attachment step always finds a sequence for each leaf.
  let sequences = complete_alignment_for_leaves(&graph, sequences, &alphabet, params.ignore_missing_alns, names)?;

  let alignment_length = get_common_length(&sequences)?;
  let mask = create_mask(&sequences, alignment_length, &alphabet);
  let mut rng = get_random_number_generator(params.seed);

  match params.method {
    MethodAncestral::Parsimony => {
      progress.check_cancelled()?;
      progress.report("Fitch parsimony", 0.3, "");
      let partition = create_fitch_partition(&graph, 0, alphabet, &sequences, names)?;
      let mut partitions_parsimony = vec![partition];

      if params.impute_missing_data {
        log::warn!(
          "--impute-missing-data has no effect with --method-anc=parsimony: Fitch parsimony produces no \
           posterior profile to impute missing tip states from. Leaf states are emitted as observed."
        );
      }

      let node_sequences =
        ancestral_reconstruction_fitch(&graph, params.include_leaves, &mut partitions_parsimony, |node, seq| {
          on_sequence(node.key, seq)
        })?;

      let partition = partitions_parsimony
        .into_iter()
        .next()
        .expect("partition vec not empty");
      progress.report("Done", 1.0, "");
      Ok(AncestralOutputFull {
        output: AncestralOutput {
          graph,
          gtr: None,
          model_name: params.model,
          mask,
          node_sequences,
        },
        partition: Some(AncestralPartition::Fitch(partition)),
      })
    },
    MethodAncestral::Marginal => {
      progress.check_cancelled()?;
      progress.report("Inferring GTR model", 0.2, "");

      let created = create_marginal_partition(
        &graph,
        0,
        alphabet,
        &sequences,
        params.model,
        params.dense,
        branch_lengths,
        names,
      )?;
      let model_name = created.model_name;
      let refine = params.gtr_iterations > 0 && params.model == GtrModelName::Infer;

      match created.partition {
        MarginalPartition::Sparse(partition, node_states) => {
          progress.check_cancelled()?;
          progress.report("Marginal reconstruction", 0.4, "");
          let (node_states, backward, forward, estimates, _log_lh) =
            partition.marginal_update(&graph, &profile_lengths, node_states)?;

          let (partition, mut node_states, backward, forward, estimates) = if refine {
            refine_gtr_iterative(
              &graph,
              partition,
              branch_lengths,
              node_states,
              backward,
              forward,
              params.gtr_iterations,
              None,
              1.0,
              None,
              false,
            )
            .map(|(p, n, b, f, e, _lh)| (p, n, b, f, e))?
          } else {
            (partition, node_states, backward, forward, estimates)
          };

          progress.check_cancelled()?;
          progress.report("Reconstructing sequences", 0.6, "");
          let node_sequences = ancestral_reconstruction(
            &graph,
            |node| {
              partition.reconstruct_node_sequence(
                &mut node_states,
                &forward,
                node,
                params.include_leaves,
                params.impute_missing_data,
                params.sample_from_profile,
                &mut rng,
              )
            },
            |key, seq| on_sequence(key, seq),
          )?;

          let gtr = partition.gtr().clone();
          progress.report("Done", 1.0, "");
          Ok(AncestralOutputFull {
            output: AncestralOutput {
              graph,
              gtr: Some(gtr),
              model_name,
              mask,
              node_sequences,
            },
            partition: Some(AncestralPartition::Sparse(SparseReconstruction {
              partition,
              node_states,
              backward,
              forward,
              estimates,
            })),
          })
        },
        MarginalPartition::Dense(partition) => {
          progress.check_cancelled()?;
          progress.report("Marginal reconstruction", 0.4, "");
          let node_states = partition.attach_sequences(&graph, &sequences, names)?;
          let (node_states, backward, forward, estimates, _log_lh) =
            partition.marginal_update(&graph, &profile_lengths, node_states)?;

          let (partition, mut node_states, backward, forward, estimates) = if refine {
            refine_gtr_iterative(
              &graph,
              partition,
              branch_lengths,
              node_states,
              backward,
              forward,
              params.gtr_iterations,
              None,
              1.0,
              None,
              false,
            )
            .map(|(p, n, b, f, e, _lh)| (p, n, b, f, e))?
          } else {
            (partition, node_states, backward, forward, estimates)
          };

          progress.check_cancelled()?;
          progress.report("Reconstructing sequences", 0.6, "");
          let node_sequences = ancestral_reconstruction(
            &graph,
            |node| {
              partition.reconstruct_node_sequence(
                &mut node_states,
                node,
                params.include_leaves,
                params.impute_missing_data,
                params.sample_from_profile,
                &mut rng,
              )
            },
            |key, seq| on_sequence(key, seq),
          )?;

          let gtr = partition.gtr().clone();
          progress.report("Done", 1.0, "");
          Ok(AncestralOutputFull {
            output: AncestralOutput {
              graph,
              gtr: Some(gtr),
              model_name,
              mask,
              node_sequences,
            },
            partition: Some(AncestralPartition::Dense(DenseReconstruction {
              partition,
              node_states,
              backward,
              forward,
              estimates,
            })),
          })
        },
      }
    },
    MethodAncestral::Joint => {
      let available = MethodAncestral::VARIANTS
        .iter()
        .filter(|v| **v != "joint")
        .copied()
        .collect::<Vec<_>>()
        .join(", ");
      make_error!("Joint ancestral reconstruction has been removed. Available methods: {available}")
    },
  }
}
