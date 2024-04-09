use crate::commands::ancestral::anc_args::TreetimeAncestralArgs;
use crate::commands::ancestral::anc_graph::{create_graph, infer_graph};
use crate::io::fasta::{read_many_fasta, FastaRecord, FastaWriter};
use crate::io::file::create_file;
use crate::io::nex::{write_nex, WriteNexOptions};
use crate::io::nwk::{write_nwk_writer, WriteNwkOptions};
use crate::seq::representation::{compress_sequences, reconstruct_ancestral_sequences};
use crate::utils::random::get_random_number_generator;
use crate::utils::string::vec_to_string;
use eyre::Report;

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

  let graph = match tree {
    None => infer_graph()?,
    Some(tree) => create_graph(tree)?,
  };

  graph.print_graph(create_file(outdir.join("graph_input.dot"))?)?;

  // TODO: avoid reading all sequences into memory somehow?
  let fastas = read_many_fasta(input_fastas)?;
  let seqs = fastas
    .into_iter()
    .map(|FastaRecord { seq_name, seq, .. }| (seq_name, seq))
    .collect();

  compress_sequences(&seqs, &graph, &mut rng)?;

  let fasta_file = create_file(outdir.join("ancestral_sequences.fasta"))?;
  let mut fasta_writer = FastaWriter::new(fasta_file);
  reconstruct_ancestral_sequences(&graph, *reconstruct_tip_states, |node, seq| {
    // TODO: avoid converting vec to string, write vec chars directly
    fasta_writer.write(&node.name, &vec_to_string(seq.to_owned())).unwrap();
  })?;

  graph.print_graph(create_file(outdir.join("graph_output.dot"))?)?;

  write_nwk_writer(
    &mut create_file(outdir.join("annotated_tree.nwk"))?,
    &graph,
    &WriteNwkOptions::default(),
  )?;

  write_nex(
    &mut create_file(outdir.join("annotated_tree.nexus"))?,
    &graph,
    &WriteNexOptions::default(),
  )?;

  Ok(())
}
