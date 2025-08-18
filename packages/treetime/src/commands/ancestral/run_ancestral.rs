use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::anc_args::{MethodAncestral, TreetimeAncestralArgs};
use crate::commands::ancestral::fitch::{ancestral_reconstruction_fitch, compress_sequences};
use crate::commands::ancestral::marginal_dense::{ancestral_reconstruction_marginal_dense, run_marginal_dense};
use crate::commands::ancestral::marginal_sparse::{ancestral_reconstruction_marginal_sparse, run_marginal_sparse};
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{Described, GraphNode, Named};
use crate::gtr::get_gtr::{get_gtr, get_gtr_dense};
use crate::io::fasta::{FastaReader, FastaRecord, FastaWriter, read_many_fasta};
use crate::io::file::{create_file_or_stdout, open_stdin};
use crate::io::nex::{NexWriteOptions, nex_write_file};
use crate::io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use crate::representation::infer_dense::infer_dense;
use crate::representation::partitions_likelihood::{PartitionLikelihood, PartitionLikelihoodWithAln};
use crate::representation::partitions_parsimony::{PartitionParsimony, PartitionParsimonyWithAln};
use crate::representation::raw_partition::RawPartition;
use crate::representation::repr_graph::{ReprGraph, ReprNode};
use crate::representation::seq::Seq;
use crate::utils::random::get_random_number_generator;
use eyre::Report;
use itertools::Itertools;
use log::info;
use std::path::Path;

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
    model_name,
    gtr_params,
    aa,
    keep_overhangs,
    zero_based,
    reconstruct_tip_states,
    report_ambiguous,
    method_anc,
    dense,
    outdir,
    seed,
  } = ancestral_args;

  let rng = get_random_number_generator(*seed);

  let dense = dense.unwrap_or_else(infer_dense);

  let treat_gap_as_unknown = dense;
  let alphabet = Alphabet::new(alphabet.unwrap_or_default(), treat_gap_as_unknown)?;

  // TODO: avoid reading all sequences into memory somehow?
  let aln = if input_fastas.is_empty() {
    info!("Reading input fasta from standard input");
    let mut reader = FastaReader::new(open_stdin()?, &alphabet);
    let mut record = FastaRecord::default();
    reader.read(&mut record)?;
    vec![record]
  } else {
    read_many_fasta(input_fastas, &alphabet)?
  };

  let output_fasta = create_file_or_stdout(outdir.join("ancestral_sequences.fasta"))?;
  let mut output_fasta = FastaWriter::new(output_fasta);

  let graph: ReprGraph = nwk_read_file(tree)?;

  let partitions = match method_anc {
    MethodAncestral::Parsimony => {
      vec![RawPartition::parsimony(alphabet.clone(), aln)?]
    },
    MethodAncestral::Marginal => {
      if !dense {
        let gtr = get_gtr(model_name, &alphabet, &graph)?;
        vec![RawPartition::likelihood(gtr, alphabet.clone(), aln)?]
      } else {
        let gtr = get_gtr_dense(model_name, &alphabet, &graph)?;
        vec![RawPartition::likelihood(gtr, alphabet.clone(), aln)?]
      }
    },
    MethodAncestral::Joint => {
      unimplemented!("MethodAncestral::MaximumLikelihoodJoint")
    },
  };

  run_inference(&graph, &partitions);

  run_reconstruction(&graph, *reconstruct_tip_states, &partitions, |node, seq| {
    let name = node.get_name_maybe().unwrap_or_default();
    let desc = node.get_desc();
    output_fasta.write(name, &desc, seq).unwrap();
  })?;

  write_graph(outdir, &graph)?;

  match method_anc {
    // both MaximumLikelihoodJoint and MaximumLikelihoodMarginal need an GTR, parsimony does not
    // it thus might make sense to split parsimony from the other two methods. For parsimony, we
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
    // only sensible for the probabilistic methods. We can deal with this later, but good to keep in mind.

    // VCF input is basically another way to instantiate the sparse representation. MAT format would be another one.
    // Again, we can deal with this later, but the entrypoint is not always going to be a fasta file.
    MethodAncestral::Parsimony => {
      let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
      let partitions = compress_sequences(&graph, partitions)?;

      ancestral_reconstruction_fitch(&graph, *reconstruct_tip_states, &partitions, |node, seq| {
        let name = node.get_name_maybe().unwrap_or_default();
        let desc = node.get_desc();
        output_fasta.write(name, &desc, seq).unwrap();
      })?;
    },
    MethodAncestral::Marginal => {
      if !dense {
        let partitions = vec![PartitionParsimonyWithAln::new(alphabet.clone(), aln)?];
        let partitions = compress_sequences(&graph, partitions)?;

        let gtr = get_gtr(model_name, &alphabet, &graph)?;
        let partitions = partitions
          .into_iter()
          .map(|part| PartitionLikelihood::from_parsimony(gtr.clone(), part)) // FIXME: avoid cloning
          .collect_vec();

        run_marginal_sparse(&graph, &partitions)?;

        ancestral_reconstruction_marginal_sparse(&graph, *reconstruct_tip_states, &partitions, |node, seq| {
          let name = node.get_name_maybe().unwrap_or_default();
          let desc = node.get_desc();
          output_fasta.write(name, &desc, &seq).unwrap();
        })?;
      } else {
        let gtr = get_gtr_dense(model_name, &alphabet, &graph)?;

        let partitions = vec![PartitionLikelihoodWithAln::new(gtr, alphabet, aln)?];
        run_marginal_dense(&graph, partitions, true)?;

        ancestral_reconstruction_marginal_dense(&graph, *reconstruct_tip_states, |node, seq| {
          let name = node.get_name().unwrap_or_default();
          let desc = node.get_desc();
          output_fasta.write(name, &desc, seq).unwrap();
        })?;
      }
    },
    MethodAncestral::Joint => {
      unimplemented!("MethodAncestral::MaximumLikelihoodJoint")
    },
  }

  Ok(())
}

pub fn run_inference(graph: &ReprGraph, partitions: &[RawPartition]) {}

pub fn run_reconstruction(
  graph: &ReprGraph,
  include_leaves: bool,
  partitions: &[RawPartition],
  mut visitor: impl FnMut(&ReprNode, &Seq),
) -> Result<(), Report> {
  Ok(())
}

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Send + Sync + Default,
{
  // json_write_file(
  //   outdir.as_ref().join("annotated_tree.graph.json"),
  //   &graph,
  //   JsonPretty(true),
  // )?;

  nwk_write_file(
    outdir.as_ref().join("annotated_tree.nwk"),
    graph,
    &NwkWriteOptions::default(),
  )?;

  nex_write_file(
    outdir.as_ref().join("annotated_tree.nexus"),
    graph,
    &NexWriteOptions::default(),
  )?;

  Ok(())
}
