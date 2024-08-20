use crate::commands::ancestral::anc_args::{MethodAncestral, TreetimeAncestralArgs};
use crate::commands::ancestral::anc_graph::{create_graph, infer_graph};
use crate::commands::ancestral::anc_reconstruction_fitch::ancestral_reconstruction_fitch;
use crate::commands::ancestral::anc_reconstruction_marginal::ancestral_reconstruction_marginal;
use crate::gtr::get_gtr::get_gtr;
use crate::io::fasta::{read_many_fasta, FastaRecord, FastaWriter};
use crate::io::file::create_file_or_stdout;
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
  json_write_file(outdir.join("graph_input.json"), &graph, JsonPretty(true))?;

  // TODO: avoid reading all sequences into memory somehow?
  let fastas = read_many_fasta(input_fastas)?;
  let seqs = fastas
    .into_iter()
    .map(|FastaRecord { seq_name, seq, .. }| (seq_name, seq))
    .collect();

  // we might want to include a heuristic when to use the dense or sparse representation
  // generally, long branches in the tree --> dense, short branches --> sparse
  // for small datasets, it doesn't really matter, but for large ones the sparse is more memory efficient when branches are short
  compress_sequences(&seqs, &graph, &mut rng)?;

  let fasta_file = create_file_or_stdout(outdir.join("ancestral_sequences.fasta"))?;
  let mut fasta_writer = FastaWriter::new(fasta_file);

  match method_anc {
    // both MaximumLikelihoodJoint and MaximumLikelihoodMarginal need an GTR, parsimony does not
    // it thus might make sense to split parsimony from the othe two methods. For parsimony, we
    // always compress the sequences and we are done (sparse in the new code). For the other two
    // there is the split between in dense and sparse.

    // To infer a GTR from the data, we can use the compressed sequence representation
    // and calculate the total number of mutations for each pair of characters, as well as the
    // time spent on each character. This information can be obtained by iterating over edges
    // to count mutations from j to i as n[i,j] and the time spent in i as t[i].
    // To calculate the t[i] for each nucleotide i, we multiply the branch_length and the count of
    // the nucleotide in the sequence composition "fixed_counts". The root_state is the same as the fixed_counts of the root.

    // One thing we haven't dealt with at all yet is how to create partitions. One common use case
    // I imagine is partition by codon position. This would require an annotation (gff). Again, this is
    // only sensible for the probabilitistic methods. We can deal with this later, but good to keep in mind.

    // VCF input is basically another way to instantiate the sparse representation. MAT format would be another one.
    // Again, we can deal with this later, but the entrypoint is not always going to be a fasta file.
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
        let name = node.name.as_deref().unwrap_or("");
        // TODO: avoid converting vec to string, write vec chars directly
        fasta_writer.write(name, &vec_to_string(seq.to_owned())).unwrap();
      })?;
    }
    MethodAncestral::Parsimony => {
      ancestral_reconstruction_fitch(&graph, *reconstruct_tip_states, |node, seq| {
        let name = node.name.as_deref().unwrap_or("");
        // TODO: avoid converting vec to string, write vec chars directly
        fasta_writer.write(name, &vec_to_string(seq.to_owned())).unwrap();
      })?;
    }
  }

  graphviz_write_file(outdir.join("graph_output.dot"), &graph)?;
  json_write_file(outdir.join("graph_output.json"), &graph, JsonPretty(true))?;

  nwk_write(
    &mut create_file_or_stdout(outdir.join("annotated_tree.nwk"))?,
    &graph,
    &NwkWriteOptions::default(),
  )?;

  nex_write(
    &mut create_file_or_stdout(outdir.join("annotated_tree.nexus"))?,
    &graph,
    &NexWriteOptions::default(),
  )?;

  Ok(())
}
