use crate::alphabet::alphabet::Alphabet;
use crate::alphabet::find_mutations::find_mutations;
use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::anc_args::TreetimeAncestralArgs;
use crate::ancestral::anc_graph::{create_graph, infer_graph, AncestralGraph};
use crate::ancestral::anc_seq::reconstruct_ancestral_sequences;
use crate::ancestral::run_anc_method::run_anc_method;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::{GraphNodeForward, NodeEdgePayloadPair};
use crate::gtr::get_gtr::get_gtr;
use crate::io::fasta::{FastaRecord, FastaWriter};
use crate::io::file::create_file;
use crate::io::nex::write_nex;
use crate::make_internal_error;
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

  let mut sequence_data = SequenceData::new(input_fastas, alphabet.ambiguous())?;

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

  let fasta_records = reconstruct_ancestral_sequences(&mut sequence_data, &graph, ancestral_args, &ancestral_params);

  find_graph_mutations(&mut graph, &sequence_data);

  let fasta_file = create_file(outdir.join("ancestral_sequences.fasta"))?;
  let mut fasta_writer = FastaWriter::new(fasta_file);

  for FastaRecord { seq, seq_name, .. } in fasta_records {
    fasta_writer.write(&seq_name, &seq)?;
  }

  graph.print_graph(create_file(outdir.join("graph_output.dot"))?)?;

  let model_string = model.to_string();
  let mut model_file = create_file(outdir.join("sequence_evolution_model.txt"))?;
  writeln!(model_file, "{model_string}")?;
  info!("{model_string}");

  write_nex(&mut create_file(outdir.join("annotated_tree.nexus"))?, &graph)?;

  Ok(())
}

/// Finds mutations on a graph and attach them to node payloads.
fn find_graph_mutations(graph: &mut AncestralGraph, sequence_data: &SequenceData) {
  graph.par_iter_breadth_first_forward(
    |GraphNodeForward {
       is_root,
       is_leaf,
       key,
       payload: node,
       parents,
     }| {
      if is_root {
        return GraphTraversalContinuation::Continue;
      }

      if parents.len() > 1 {
        unimplemented!("Multiple parent nodes are not supported yet");
      }

      let (parent, _) = &parents[0];

      let parent_seq = sequence_data.get_full(&parent.read().name).unwrap();
      let this_seq = sequence_data.get_full(&node.name).unwrap();
      node.mutations = find_mutations(&parent_seq, &this_seq, &node.mask);

      GraphTraversalContinuation::Continue
    },
  );
}
