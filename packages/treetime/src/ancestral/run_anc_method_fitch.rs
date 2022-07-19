use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::anc_graph::AncestralGraph;
use crate::ancestral::run_ancestral_reconstruction::TreetimeAncestralParams;
use crate::cli::treetime_cli::TreetimeAncestralArgs;
use crate::gtr::gtr::GTR;
use eyre::Report;
use rand::Rng;

pub fn run_anc_method_fitch(
  sequence_data: &SequenceData,
  alphabet: &Alphabet,
  model: &GTR,
  graph: &mut AncestralGraph,
  rng: &mut impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  unimplemented!("ml_anc_fitch: not yet implemented");
}
