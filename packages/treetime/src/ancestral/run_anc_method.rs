use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::anc_graph::AncestralGraph;
use crate::ancestral::run_anc_method_fitch::run_anc_method_fitch;
use crate::ancestral::run_anc_method_ml_joint::run_anc_method_ml_joint;
use crate::ancestral::run_anc_method_ml_marginal::run_anc_method_ml_marginal;
use crate::ancestral::run_ancestral_reconstruction::TreetimeAncestralParams;
use crate::cli::treetime_cli::{MethodAncestral, TreetimeAncestralArgs};
use crate::gtr::gtr::GTR;
use eyre::Report;
use rand::Rng;

pub fn run_anc_method(
  sequence_data: &SequenceData,
  alphabet: &Alphabet,
  model: &GTR,
  graph: &mut AncestralGraph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  match ancestral_args.method_anc {
    MethodAncestral::MaximumLikelihoodJoint => run_anc_method_ml_joint(
      sequence_data,
      alphabet,
      model,
      graph,
      rng,
      ancestral_args,
      ancestral_params,
    ),
    MethodAncestral::MaximumLikelihoodMarginal => run_anc_method_ml_marginal(
      sequence_data,
      alphabet,
      model,
      graph,
      rng,
      ancestral_args,
      ancestral_params,
    ),
    MethodAncestral::Parsimony => run_anc_method_fitch(
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
