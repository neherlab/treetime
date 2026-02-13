use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::args::{MethodAncestral, TreetimeAncestralArgs};
use crate::commands::ancestral::fitch::{ancestral_reconstruction_fitch, compress_sequences, get_common_length};
use crate::commands::ancestral::marginal::{ancestral_reconstruction_marginal, initialize_marginal, update_marginal};
use crate::gtr::get_gtr::{JC69Params, get_gtr, jc69};
use treetime_io::fasta::{FastaReader, FastaRecord, FastaWriter, read_many_fasta};
use treetime_io::nex::{NexWriteOptions, nex_write_file};
use treetime_io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::infer_dense::infer_dense;
use crate::representation::partition_fitch::PartitionFitch;
use crate::representation::partition_marginal_dense::PartitionMarginalDense;
use crate::representation::partition_marginal_sparse::PartitionMarginalSparse;
use eyre::Report;
use itertools::Itertools;
use log::info;
use maplit::btreemap;
use parking_lot::RwLock;
use serde::Serialize;
use std::path::Path;
use std::sync::Arc;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::io::file::{create_file_or_stdout, open_stdin};

#[derive(Clone, Debug, Default)]
pub struct TreetimeAncestralParams {
  pub sample_from_profile: bool,
  pub fixed_pi: bool,
}

pub fn run_ancestral_reconstruction(ancestral_args: &TreetimeAncestralArgs) -> Result<(), Report> {
  let TreetimeAncestralArgs {
    input_fastas,
    tree,
    alphabet,
    model_name,
    reconstruct_tip_states,
    method_anc,
    dense,
    outdir,
    ..
  } = ancestral_args;

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

  let graph: GraphAncestral = nwk_read_file(tree)?;

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
      #[allow(clippy::iter_on_single_items)]
      let partitions_parsimony = [PartitionFitch {
        index: 0,
        alphabet,
        length: get_common_length(&aln)?,
        nodes: btreemap! {},
        edges: btreemap! {},
      }]
      .into_iter()
      .map(|p| Arc::new(RwLock::new(p)))
      .collect_vec();

      if !partitions_parsimony.is_empty() {
        compress_sequences(&graph, &partitions_parsimony, &aln)?;

        ancestral_reconstruction_fitch(&graph, *reconstruct_tip_states, &partitions_parsimony, |node, seq| {
          let name = node.payload.name.as_deref().unwrap_or("");
          let desc = &node.payload.desc;
          output_fasta.write(name, desc, seq).unwrap();
        })?;
      }
    },
    MethodAncestral::Marginal => {
      if !dense {
        #[allow(clippy::iter_on_single_items)]
        let partitions_marginal_sparse = [PartitionMarginalSparse {
          index: 0,
          gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
          alphabet,
          length: get_common_length(&aln)?,
          nodes: btreemap! {},
          edges: btreemap! {},
        }]
        .into_iter()
        .map(|p| Arc::new(RwLock::new(p)))
        .collect_vec();

        if !partitions_marginal_sparse.is_empty() {
          compress_sequences(&graph, &partitions_marginal_sparse, &aln)?;

          // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
          // FIXME: spaghetti code: dummy gtr is replaced by real gtr here
          for partition in &partitions_marginal_sparse {
            let gtr = get_gtr(model_name, partition, &graph)?;
            partition.write_arc().gtr = gtr;
          }

          update_marginal(&graph, &partitions_marginal_sparse)?;
          ancestral_reconstruction_marginal(
            &graph,
            *reconstruct_tip_states,
            &partitions_marginal_sparse,
            |node, seq| {
              let name = node.name.as_deref().unwrap_or("");
              let desc = &node.desc;
              output_fasta.write(name, desc, seq).unwrap();
            },
          )?;
        }
      } else {
        #[allow(clippy::iter_on_single_items)]
        let partitions_marginal_dense = [PartitionMarginalDense {
          index: 0,
          gtr: jc69(JC69Params::default())?, // FIXME: dummy temporary gtr should not be needed here
          alphabet,
          length: get_common_length(&aln)?,
          nodes: btreemap! {},
          edges: btreemap! {},
        }]
        .into_iter()
        .map(|p| Arc::new(RwLock::new(p)))
        .collect_vec();

        if !partitions_marginal_dense.is_empty() {
          initialize_marginal(&graph, &partitions_marginal_dense, &aln)?;

          // FIXME: chicken & egg problem: to get a gtr we need partitions, to get partitions we need a gtr
          // FIXME: spaghetti code: dummy gtr is replaced by real gtr here
          for partition in &partitions_marginal_dense {
            let gtr = get_gtr(model_name, partition, &graph)?;
            partition.write_arc().gtr = gtr;
          }

          ancestral_reconstruction_marginal(
            &graph,
            *reconstruct_tip_states,
            &partitions_marginal_dense,
            |node, seq| {
              let name = node.name.as_deref().unwrap_or("");
              let desc = &node.desc;
              output_fasta.write(name, desc, seq).unwrap();
            },
          )?;
        }
      }
    },
    MethodAncestral::Joint => {
      unimplemented!("MethodAncestral::MaximumLikelihoodJoint")
    },
  }

  write_graph(outdir, &graph)?;

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
