use crate::commands::ancestral::anc_args::{MethodAncestral, TreetimeAncestralArgs};
use crate::commands::ancestral::anc_graph::{create_graph, infer_graph};
use crate::commands::ancestral::anc_reconstruction_fitch::ancestral_reconstruction_fitch;
use crate::commands::ancestral::anc_reconstruction_marginal::ancestral_reconstruction_marginal;
use crate::gtr::get_gtr::get_gtr;
use crate::io::fasta::{read_many_fasta, FastaRecord, FastaWriter};
use crate::io::file::create_file;
use crate::io::graphviz::graphviz_write_file;
use crate::io::json::{json_write_file, JsonPretty};
use crate::io::nex::{nex_write, NexWriteOptions};
use crate::io::nwk::{nwk_write, NwkWriteOptions};
use crate::seq::representation::compress_sequences;
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

  graphviz_write_file(outdir.join("graph_input.dot"), &graph)?;
  json_write_file(&graph, outdir.join("graph_input.json"), JsonPretty(true))?;

  // TODO: avoid reading all sequences into memory somehow?
  let fastas = read_many_fasta(input_fastas)?;
  let seqs = fastas
    .into_iter()
    .map(|FastaRecord { seq_name, seq, .. }| (seq_name, seq))
    .collect();

  compress_sequences(&seqs, &graph, &mut rng)?;

  let fasta_file = create_file(outdir.join("ancestral_sequences.fasta"))?;
  let mut fasta_writer = FastaWriter::new(fasta_file);

  match method_anc {
    MethodAncestral::MaximumLikelihoodJoint => {
      unimplemented!("MethodAncestral::MaximumLikelihoodJoint")
    }
    MethodAncestral::MaximumLikelihoodMarginal => {
      let gtr = get_gtr(gtr, alphabet)?;

      // Uncomment this for custom GTR
      //
      // let alphabet = Alphabet::new(AlphabetName::NucNogap)?;
      //
      // let gtr = GTR::new(GTRParams {
      //   alphabet,
      //   mu: 1.0,
      //   W: None,
      //   pi: array![0.2, 0.3, 0.15, 0.45],
      // })?;

      ancestral_reconstruction_marginal(&graph, &gtr, *reconstruct_tip_states, |node, seq| {
        // TODO: avoid converting vec to string, write vec chars directly
        fasta_writer.write(&node.name, &vec_to_string(seq.to_owned())).unwrap();
      })?;
    }
    MethodAncestral::Parsimony => {
      ancestral_reconstruction_fitch(&graph, *reconstruct_tip_states, |node, seq| {
        // TODO: avoid converting vec to string, write vec chars directly
        fasta_writer.write(&node.name, &vec_to_string(seq.to_owned())).unwrap();
      })?;
    }
  }

  graphviz_write_file(outdir.join("graph_output.dot"), &graph)?;
  json_write_file(&graph, outdir.join("graph_output.json"), JsonPretty(true))?;

  nwk_write(
    &mut create_file(outdir.join("annotated_tree.nwk"))?,
    &graph,
    &NwkWriteOptions::default(),
  )?;

  nex_write(
    &mut create_file(outdir.join("annotated_tree.nexus"))?,
    &graph,
    &NexWriteOptions::default(),
  )?;

  Ok(())
}
