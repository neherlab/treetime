use crate::alphabet::alphabet::Alphabet;
use crate::ancestral::fitch::{ancestral_reconstruction_fitch, compress_sequences};
use crate::ancestral::marginal::{ancestral_reconstruction_marginal, initialize_marginal, update_marginal};
use crate::ancestral::mask::create_mask;
use crate::ancestral::params::MethodAncestral;
use crate::ancestral::sample::SampleMode;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::gtr::refinement::refine_gtr_iterative;
use crate::partition::create::{MarginalPartition, create_marginal_partition};
use crate::partition::fitch::PartitionFitch;
use crate::partition::marginal_dense::PartitionMarginalDense;
use crate::partition::traits::HasGtr;
use crate::payload::ancestral::{GraphAncestral, NodeAncestral};
use crate::progress::ProgressSink;
use crate::seq::alignment::get_common_length;
use eyre::Report;
use maplit::btreemap;
use parking_lot::RwLock;
use serde::Serialize;
use std::sync::Arc;
use strum::VariantNames;
use treetime_io::fasta::FastaRecord;
use treetime_primitives::Seq;
use treetime_utils::make_error;
use treetime_utils::sync::random::get_random_number_generator;

pub struct AncestralParams {
  pub method: MethodAncestral,
  pub model: GtrModelName,
  pub dense: Option<bool>,
  pub reconstruct_tip_states: bool,
  pub gtr_iterations: usize,
  pub site_specific_gtr: bool,
  pub seed: Option<u64>,
  pub sample_from_profile: SampleMode,
}

pub struct AncestralInput {
  pub graph: GraphAncestral,
  pub alphabet: Alphabet,
  pub sequences: Vec<FastaRecord>,
}

/// Partition variant returned from ancestral reconstruction, allowing the caller to build comment providers.
pub enum AncestralPartition {
  Fitch(Arc<RwLock<PartitionFitch>>),
  Sparse(Arc<RwLock<crate::partition::marginal_sparse::PartitionMarginalSparse>>),
  Dense(Arc<RwLock<PartitionMarginalDense>>),
}

#[derive(Debug, Serialize)]
pub struct AncestralOutput {
  #[serde(skip)]
  pub graph: GraphAncestral,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub mask: Vec<bool>,
}

pub struct AncestralOutputFull {
  pub output: AncestralOutput,
  pub partition: Option<AncestralPartition>,
}

pub fn run<F>(
  params: &AncestralParams,
  input: AncestralInput,
  mut on_sequence: F,
  progress: &dyn ProgressSink,
) -> Result<AncestralOutputFull, Report>
where
  F: FnMut(&NodeAncestral, &Seq) -> Result<(), Report>,
{
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

  let alignment_length = get_common_length(&input.sequences)?;
  let mask = create_mask(&input.sequences, alignment_length, &input.alphabet);
  let mut rng = get_random_number_generator(params.seed);

  match params.method {
    MethodAncestral::Parsimony => {
      progress.check_cancelled()?;
      progress.report("Fitch parsimony", 0.3, "");
      let partition = Arc::new(RwLock::new(PartitionFitch {
        index: 0,
        alphabet: input.alphabet,
        length: alignment_length,
        nodes: btreemap! {},
        edges: btreemap! {},
      }));
      let partitions_parsimony = vec![Arc::clone(&partition)];

      compress_sequences(&input.graph, &partitions_parsimony, &input.sequences)?;

      ancestral_reconstruction_fitch(&input.graph, params.reconstruct_tip_states, &partitions_parsimony, |node, seq| {
        on_sequence(&node.payload, seq)
      })?;

      progress.report("Done", 1.0, "");
      Ok(AncestralOutputFull {
        output: AncestralOutput {
          graph: input.graph,
          gtr: None,
          model_name: params.model,
          mask,
        },
        partition: Some(AncestralPartition::Fitch(partition)),
      })
    },
    MethodAncestral::Marginal => {
      progress.check_cancelled()?;
      progress.report("Inferring GTR model", 0.2, "");

      let created = create_marginal_partition(&input.graph, 0, input.alphabet, &input.sequences, params.model, params.dense)?;

      match created.partition {
        MarginalPartition::Sparse(partition) => {
          let partitions = vec![Arc::new(RwLock::new(partition))];

          progress.check_cancelled()?;
          progress.report("Marginal reconstruction", 0.4, "");
          update_marginal(&input.graph, &partitions)?;

          if params.gtr_iterations > 0 && params.model == GtrModelName::Infer {
            refine_gtr_iterative(&input.graph, &partitions[0], params.gtr_iterations, None, 1.0, None, false)?;
          }

          progress.check_cancelled()?;
          progress.report("Reconstructing sequences", 0.6, "");
          ancestral_reconstruction_marginal(
            &input.graph,
            params.reconstruct_tip_states,
            &partitions,
            params.sample_from_profile,
            &mut rng,
            |node, seq| on_sequence(node, seq),
          )?;

          let gtr = partitions[0].read_arc().gtr().clone();
          progress.report("Done", 1.0, "");
          Ok(AncestralOutputFull {
            output: AncestralOutput {
              graph: input.graph,
              gtr: Some(gtr),
              model_name: created.model_name,
              mask,
            },
            partition: Some(AncestralPartition::Sparse(partitions.into_iter().next().expect("partition vec not empty"))),
          })
        },
        MarginalPartition::Dense(partition) => {
          let partitions = vec![Arc::new(RwLock::new(partition))];

          progress.check_cancelled()?;
          progress.report("Marginal reconstruction", 0.4, "");
          initialize_marginal(&input.graph, &partitions, &input.sequences)?;
          update_marginal(&input.graph, &partitions)?;

          if params.gtr_iterations > 0 && params.model == GtrModelName::Infer {
            refine_gtr_iterative(&input.graph, &partitions[0], params.gtr_iterations, None, 1.0, None, false)?;
          }

          progress.check_cancelled()?;
          progress.report("Reconstructing sequences", 0.6, "");
          ancestral_reconstruction_marginal(
            &input.graph,
            params.reconstruct_tip_states,
            &partitions,
            params.sample_from_profile,
            &mut rng,
            |node, seq| on_sequence(node, seq),
          )?;

          let gtr = partitions[0].read_arc().gtr().clone();
          progress.report("Done", 1.0, "");
          Ok(AncestralOutputFull {
            output: AncestralOutput {
              graph: input.graph,
              gtr: Some(gtr),
              model_name: created.model_name,
              mask,
            },
            partition: Some(AncestralPartition::Dense(partitions.into_iter().next().expect("partition vec not empty"))),
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
