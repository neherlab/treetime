use crate::alphabet::sequence_data::SequenceData;
use crate::commands::ancestral::anc_args::TreetimeAncestralArgs;
use crate::commands::ancestral::anc_graph::AncestralGraph;
use crate::commands::ancestral::run_ancestral_reconstruction::TreetimeAncestralParams;
use crate::gtr::gtr::GTR;
use eyre::Report;
use rand::Rng;

pub fn run_anc_method_ml_marginal(
  sequence_data: &SequenceData,
  model: &GTR,
  graph: &AncestralGraph,
  rng: &impl Rng,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Result<(), Report> {
  unimplemented!("ml_anc_marginal: not yet implemented");
}