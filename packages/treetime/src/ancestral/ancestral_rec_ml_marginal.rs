use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::ancestral_graph::Graph;
use crate::ancestral::run_ancestral::TreetimeAncestralParams;
use crate::cli::treetime_cli::TreetimeAncestralArgs;
use crate::gtr::gtr::GTR;
use eyre::Report;
use rand::Rng;

pub fn ml_anc_marginal(
  sequence_data: &SequenceData,
  alphabet: &Alphabet,
  model: &GTR,
  graph: &mut Graph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  unimplemented!("ml_anc_marginal: not yet implemented");
}
