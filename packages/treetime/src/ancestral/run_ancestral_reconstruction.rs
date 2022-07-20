use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::anc_graph::{create_graph, infer_graph};
use crate::ancestral::anc_seq::reconstruct_ancestral_sequences;
use crate::ancestral::run_anc_method::run_anc_method;
use crate::cli::treetime_cli::TreetimeAncestralArgs;
use crate::gtr::get_gtr::get_gtr;
use crate::io::fasta::{FastaRecord, FastaWriter};
use crate::io::file::create_file;
use crate::utils::random::get_random_number_generator;
use eyre::Report;
use itertools::Itertools;
use std::fmt::Display;

#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Debug, Default)]
pub struct TreetimeAncestralParams {
  pub sample_from_profile: bool,
  pub fixed_pi: bool,
}

pub fn run_ancestral_reconstruction(ancestral_args: &TreetimeAncestralArgs) -> Result<(), Report> {
  let TreetimeAncestralArgs {
    input_fastas,
    aln,
    vcf_reference,
    tree,
    gtr,
    gtr_params,
    aa,
    keep_overhangs,
    zero_based,
    reconstruct_tip_states,
    report_ambiguous,
    method_anc,
    outdir,
    seed,
  } = ancestral_args;

  let mut rng = get_random_number_generator(*seed);

  // TODO: alphabet is hardcoded. Make it dynamic.
  let alphabet = Alphabet::new("nuc")?;

  let sequence_data = SequenceData::new(input_fastas, alphabet.ambiguous())?;

  let model = get_gtr(gtr)?;

  let mut graph = match tree {
    None => infer_graph()?,
    Some(tree) => create_graph(tree)?,
  };

  graph.print_graph(create_file(outdir.join("graph_input.dot"))?)?;

  let ancestral_params = TreetimeAncestralParams::default();

  run_anc_method(
    &sequence_data,
    &alphabet,
    &model,
    &mut graph,
    &mut rng,
    ancestral_args,
    &ancestral_params,
  )?;

  let fasta_records = reconstruct_ancestral_sequences(&sequence_data, &graph, ancestral_args, &ancestral_params);

  let fasta_file = create_file(outdir.join("ancestral_sequences.fasta"))?;
  let mut fasta_writer = FastaWriter::new(fasta_file);

  for FastaRecord { seq, seq_name, .. } in fasta_records {
    fasta_writer.write(&seq_name, &seq)?;
  }

  graph.print_graph(create_file(outdir.join("graph_output.dot"))?)?;

  Ok(())
}