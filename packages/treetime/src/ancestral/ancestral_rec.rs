use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::ancestral_graph::Graph;
use crate::ancestral::ancestral_rec_fitch::ml_anc_fitch;
use crate::ancestral::ancestral_rec_ml_joint::ml_anc_joint;
use crate::ancestral::ancestral_rec_ml_marginal::ml_anc_marginal;
use crate::ancestral::run_ancestral::TreetimeAncestralParams;
use crate::cli::treetime_cli::{MethodAncestral, TreetimeAncestralArgs};
use crate::gtr::gtr::GTR;
use eyre::Report;
use rand::Rng;

pub fn run_ancestral_reconstruction(
  sequence_data: &SequenceData,
  alphabet: &Alphabet,
  model: &GTR,
  graph: &mut Graph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  match ancestral_args.method_anc {
    MethodAncestral::MaximumLikelihoodJoint => ml_anc_joint(
      sequence_data,
      alphabet,
      model,
      graph,
      rng,
      ancestral_args,
      ancestral_params,
    ),
    MethodAncestral::MaximumLikelihoodMarginal => ml_anc_marginal(
      sequence_data,
      alphabet,
      model,
      graph,
      rng,
      ancestral_args,
      ancestral_params,
    ),
    MethodAncestral::Parsimony => ml_anc_fitch(
      sequence_data,
      alphabet,
      model,
      graph,
      rng,
      ancestral_args,
      ancestral_params,
    ),
  }
}
