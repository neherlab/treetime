use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::anc_args::{MethodAncestral, TreetimeAncestralArgs};
use crate::commands::ancestral::fitch::{ancestral_reconstruction_fitch, compress_sequences};
use crate::commands::ancestral::marginal_dense::{ancestral_reconstruction_marginal_dense, run_marginal_dense};
use crate::commands::ancestral::marginal_sparse::{ancestral_reconstruction_marginal_sparse, run_marginal_sparse};
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use crate::gtr::get_gtr::{get_gtr, get_gtr_dense};
use crate::io::fasta::{read_many_fasta, FastaWriter};
use crate::io::file::create_file_or_stdout;
use crate::io::nex::{nex_write_file, NexWriteOptions};
use crate::io::nwk::{nwk_read_file, nwk_write_file, EdgeToNwk, NodeToNwk, NwkWriteOptions};
use crate::representation::graph_dense::DenseGraph;
use crate::representation::graph_sparse::SparseGraph;
use crate::representation::infer_dense::infer_dense;
use crate::representation::partitions_likelihood::{PartitionLikelihood, PartitionLikelihoodWithAln};
use crate::representation::partitions_parsimony::PartitionParsimonyWithAln;
use crate::utils::random::get_random_number_generator;
use crate::utils::string::vec_to_string;
use eyre::Report;
use itertools::Itertools;
use serde::Serialize;
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

  let alphabet = Alphabet::new(alphabet.unwrap_or_default(), false)?;

  // TODO: avoid reading all sequences into memory somehow?
  let aln = read_many_fasta(input_fastas, &alphabet)?;

  let output_fasta = create_file_or_stdout(outdir.join("ancestral_sequences.fasta"))?;
  let mut output_fasta = FastaWriter::new(output_fasta);

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
    MethodAncestral::Parsimony => {
      let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
      let graph: SparseGraph = nwk_read_file(tree)?;
      let partitions = compress_sequences(&graph, partitions)?;

      ancestral_reconstruction_fitch(&graph, *reconstruct_tip_states, &partitions, |node, seq| {
        let name = node.name.as_deref().unwrap_or("");
        let desc = &node.desc;
        // TODO: avoid converting vec to string, write vec chars directly
        output_fasta.write(name, desc, vec_to_string(seq)).unwrap();
      })?;

      write_graph(outdir, &graph)?;
    }
    MethodAncestral::Marginal => {
      let dense = dense.unwrap_or_else(infer_dense);

      if !dense {
        let graph: SparseGraph = nwk_read_file(tree)?;
        let partitions = vec![PartitionParsimonyWithAln::new(alphabet.clone(), aln)?];
        let partitions = compress_sequences(&graph, partitions)?;

        let gtr = get_gtr(model_name, &alphabet, &graph)?;
        let partitions = partitions
          .into_iter()
          .map(|part| PartitionLikelihood::from_parsimony(gtr.clone(), part)) // FIXME: avoid cloning
          .collect_vec();

        run_marginal_sparse(&graph, &partitions)?;

        ancestral_reconstruction_marginal_sparse(&graph, *reconstruct_tip_states, &partitions, |node, seq| {
          let name = node.name.as_deref().unwrap_or("");
          let desc = &node.desc;
          // TODO: avoid converting vec to string, write vec chars directly
          output_fasta.write(name, desc, vec_to_string(seq)).unwrap();
        })?;

        write_graph(outdir, &graph)?;
      } else {
        let graph: DenseGraph = nwk_read_file(tree)?;
        let gtr = get_gtr_dense(model_name, &alphabet, &graph)?;

        let partitions = vec![PartitionLikelihoodWithAln::new(gtr, alphabet, aln)?];
        run_marginal_dense(&graph, partitions)?;

        ancestral_reconstruction_marginal_dense(&graph, *reconstruct_tip_states, |node, seq| {
          let name = node.name.as_deref().unwrap_or("");
          let desc = &node.desc;
          // TODO: avoid converting vec to string, write vec chars directly
          output_fasta.write(name, desc, vec_to_string(seq)).unwrap();
        })?;

        write_graph(outdir, &graph)?;
      }
    }
    MethodAncestral::Joint => {
      unimplemented!("MethodAncestral::MaximumLikelihoodJoint")
    }
  };

  Ok(())
}

fn write_graph<N, E, D>(outdir: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk + Serialize,
  E: GraphEdge + EdgeToNwk + Serialize,
  D: Send + Sync + Default + Serialize,
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
