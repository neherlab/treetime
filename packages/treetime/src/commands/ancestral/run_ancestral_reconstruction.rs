use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::find_mutations::find_graph_mutations;
use crate::alphabet::sequence_data::SequenceData;
use crate::commands::ancestral::anc_args::TreetimeAncestralArgs;
use crate::commands::ancestral::anc_graph::{create_graph, infer_graph, AncestralGraph};
use crate::commands::ancestral::anc_seq::reconstruct_ancestral_sequences;
use crate::commands::ancestral::run_anc_method::run_anc_method;
use crate::gtr::get_gtr::get_gtr;
use crate::gtr::gtr::GTR;
use crate::io::fasta::{FastaRecord, FastaWriter};
use crate::io::file::create_file;
use crate::io::nex::write_nex;
use crate::utils::random::get_random_number_generator;
use eyre::Report;
use itertools::Itertools;
use log::info;
use std::fmt::Display;
use std::io::Write;

#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, Debug, Default)]
pub struct TreetimeAncestralParams {
  pub sample_from_profile: bool,
  pub fixed_pi: bool,
}

pub fn run_ancestral_reconstruction(
  ancestral_args: &TreetimeAncestralArgs,
) -> Result<(GTR, SequenceData, AncestralGraph), Report> {
  let TreetimeAncestralArgs {
    input_fastas,
    aln,
    vcf_reference,
    tree,
    alphabet,
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

  let gtr = get_gtr(gtr, alphabet)?;

  let mut sequence_data = SequenceData::new(input_fastas, gtr.alphabet().ambiguous())?;

  let mut graph = match tree {
    None => infer_graph()?,
    Some(tree) => create_graph(tree)?,
  };

  graph.print_graph(create_file(outdir.join("graph_input.dot"))?)?;

  let ancestral_params = TreetimeAncestralParams::default();

  run_anc_method(
    &sequence_data,
    &gtr,
    &mut graph,
    &mut rng,
    ancestral_args,
    &ancestral_params,
  )?;

  let fasta_records = reconstruct_ancestral_sequences(&mut sequence_data, &graph, ancestral_args, &ancestral_params);

  find_graph_mutations(&mut graph, &sequence_data);

  let fasta_file = create_file(outdir.join("ancestral_sequences.fasta"))?;
  let mut fasta_writer = FastaWriter::new(fasta_file);

  for FastaRecord { seq, seq_name, .. } in fasta_records {
    fasta_writer.write(&seq_name, &seq)?;
  }

  graph.print_graph(create_file(outdir.join("graph_output.dot"))?)?;

  let model_string = gtr.to_string();
  let mut model_file = create_file(outdir.join("sequence_evolution_model.txt"))?;
  writeln!(model_file, "{model_string}")?;
  info!("{model_string}");

  write_nex(&mut create_file(outdir.join("annotated_tree.nexus"))?, &graph)?;

  Ok((gtr, sequence_data, graph))
}
