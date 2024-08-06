use crate::alphabet::alphabet::Alphabet;
use crate::graph::node::Named;
use crate::gtr::gtr::GTR;
use crate::io::fasta::FastaRecord;
use crate::port::seq_partitions::SeqPartition;
use crate::port::seq_sparse::{SparseSeqGraph, SparseSeqNode};
use crate::{make_error, make_internal_report, make_report};
use eyre::{Report, WrapErr};
use itertools::Itertools;

#[derive(Debug)]
pub struct PartitionModel<'g> {
  aln: Vec<FastaRecord>,
  gtr: &'g GTR,
}

fn attach_seqs_to_graph<'a, 'g>(
  graph: &SparseSeqGraph<'a, 'g>,
  alphabet: &'a Alphabet,
  partitions: &[PartitionModel<'g>],
) -> Result<(), Report> {
  let mut sparse_partitions = vec![];

  for PartitionModel { aln, gtr } in partitions {
    let length = get_common_length(aln)?;
    let partition = SeqPartition::new(gtr, length, alphabet);
    sparse_partitions.push(partition);

    for leaf in graph.get_leaves() {
      let mut leaf = leaf.read_arc().payload().write_arc();

      let (leaf_name, leaf_aln) = {
        let leaf_name = leaf.name().ok_or_else(|| {
          make_report!("Expected all leaf nodes to have names, such that they can be matched with their corresponding sequences. But found a leaf node that has no name.")
        })?;
        let leaf_name = leaf_name.as_ref().to_owned();

        // TODO(perf): this might be slow if there are many sequences
        let leaf_aln = aln
          .iter()
          .find(|fasta| fasta.seq_name == leaf_name)
          .ok_or_else(|| make_internal_report!("Leaf sequence not found: '{leaf_name}'"))?;

        // TODO(perf): unnecessary copy of sequence data. Neither String, nor &[char] works well for us, it seems.
        // We probably want a custom class for sequences. Sequences should be instantiated in the fasta parser and
        // never need a copy like here.
        let leaf_aln = leaf_aln.seq.chars().collect_vec();

        (leaf_name, leaf_aln)
      };

      *leaf = SparseSeqNode::new(leaf_name, &leaf_aln, alphabet)?;
    }
  }

  graph.data().write_arc().sparse_partitions = sparse_partitions;

  Ok(())
}

pub fn get_common_length(aln: &[FastaRecord]) -> Result<usize, Report> {
  let lengths = aln
    .iter()
    .into_group_map_by(|aln| aln.seq.len())
    .into_iter()
    .collect_vec();

  match lengths[..] {
    [] => Ok(0),
    [(length, _)] => Ok(length),
    _ => {
      let message = lengths
        .into_iter()
        .sorted_by_key(|(length, _)| *length)
        .map(|(length, entries)| {
          let names = entries.iter().map(|aln| format!("    \"{}\"", aln.seq_name)).join("\n");
          format!("Length {length}:\n{names}")
        })
        .join("\n\n");

      make_error!("Sequences are expected to all have the same length, but found the following lengths:\n\n{message}")
    }
  }
  .wrap_err("When calculating length of sequences")
}
