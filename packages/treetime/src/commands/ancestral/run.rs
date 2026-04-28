use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::args::{MethodAncestral, TreetimeAncestralArgs};
use crate::commands::ancestral::fitch::{ancestral_reconstruction_fitch, compress_sequences, get_common_length};
use crate::commands::ancestral::marginal::{ancestral_reconstruction_marginal, initialize_marginal, update_marginal};
use crate::gtr::get_gtr::{JC69Params, get_gtr_dense, get_gtr_sparse, jc69, write_gtr_json};
use crate::make_error;
use crate::representation::algo::infer_dense::infer_dense;
use crate::representation::partition::fitch::PartitionFitch;
use crate::representation::partition::marginal_dense::PartitionMarginalDense;
use crate::representation::partition::marginal_sparse::PartitionMarginalSparse;
use crate::representation::partition::traits::PartitionBranchOps;
use crate::representation::payload::ancestral::{GraphAncestral, annotate_branch_mutations};
use crate::seq::gap_fill::apply_gap_fill;
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
use treetime_io::fasta::{FastaReader, FastaRecord, FastaWriter, read_many_fasta};
use treetime_io::nex::{NexWriteOptions, nex_write_file};
use treetime_io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_read_file, nwk_write_file};
use treetime_primitives::seq;
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
    site_specific_gtr,
    ..
  } = ancestral_args;

  if *site_specific_gtr {
    return make_error!(
      "--site-specific-gtr is not yet integrated into the ancestral reconstruction pipeline. \
       The mathematical core (GTRSiteSpecific) is implemented but partition system wiring is pending."
    );
  }

  let dense = dense.unwrap_or_else(infer_dense);

  let alphabet = Alphabet::new(alphabet.unwrap_or_default())?;

  let gap_fill_mode = ancestral_args.effective_gap_fill();

  // TODO: avoid reading all sequences into memory somehow?
  let mut aln = if input_fastas.is_empty() {
    info!("Reading input fasta from standard input");
    let mut reader = FastaReader::new(open_stdin()?, &alphabet);
    let mut record = FastaRecord::default();
    reader.read(&mut record)?;
    vec![record]
  } else {
    read_many_fasta(input_fastas, &alphabet)?
  };

  for record in &mut aln {
    apply_gap_fill(&mut record.seq, gap_fill_mode, alphabet.gap(), alphabet.unknown());
  }

  let output_fasta = create_file_or_stdout(outdir.join("ancestral_sequences.fasta"))?;
  let mut output_fasta = FastaWriter::new(output_fasta);

  let graph: GraphAncestral = nwk_read_file(tree)?;

  match method_anc {
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
          output_fasta.write(name, desc, seq)
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
          root_sequence: seq![],
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
            let gtr = get_gtr_sparse(model_name, partition, &graph)?;
            write_gtr_json(&gtr, *model_name, outdir, None)?;
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
              output_fasta.write(name, desc, seq)
            },
          )?;

          let branch_ops: Vec<Arc<RwLock<dyn PartitionBranchOps>>> = partitions_marginal_sparse
            .into_iter()
            .map(|p| -> Arc<RwLock<dyn PartitionBranchOps>> { p })
            .collect_vec();
          annotate_branch_mutations(&graph, &branch_ops)?;
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
          // For model inference, we need profiles populated first. Run marginal reconstruction
          // with dummy GTR, then infer the real GTR, then re-run if model changed.
          update_marginal(&graph, &partitions_marginal_dense)?;

          // Triggers second marginal pass internally when model=infer
          for partition in &partitions_marginal_dense {
            let gtr = get_gtr_dense(model_name, partition, &graph)?;
            write_gtr_json(&gtr, *model_name, outdir, None)?;
            partition.write_arc().gtr = gtr;
          }

          // Re-run with inferred GTR
          update_marginal(&graph, &partitions_marginal_dense)?;

          ancestral_reconstruction_marginal(
            &graph,
            *reconstruct_tip_states,
            &partitions_marginal_dense,
            |node, seq| {
              let name = node.name.as_deref().unwrap_or("");
              let desc = &node.desc;
              output_fasta.write(name, desc, seq)
            },
          )?;

          let branch_ops: Vec<Arc<RwLock<dyn PartitionBranchOps>>> = partitions_marginal_dense
            .into_iter()
            .map(|p| -> Arc<RwLock<dyn PartitionBranchOps>> { p })
            .collect_vec();
          annotate_branch_mutations(&graph, &branch_ops)?;
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
