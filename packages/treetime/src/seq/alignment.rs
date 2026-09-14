use crate::make_error;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_io::fasta::FastaRecord;
use treetime_io::nwk::NwkFastaNodeInput;

/// The common length of the attached leaf sequences in a node-input map, or `0` when none carry a
/// sequence. Errors when leaves disagree on length, listing each length and its leaf names.
pub fn get_common_length_of_node_inputs(
  node_inputs: &BTreeMap<GraphNodeKey, NwkFastaNodeInput>,
) -> Result<usize, Report> {
  let lengths = node_inputs
    .values()
    .filter_map(|node| node.aln.as_ref().map(|seq| (seq.len(), node)))
    .into_group_map_by(|(length, _)| *length)
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
          let names = entries
            .iter()
            .map(|(_, node)| format!("    \"{}\"", node.name.as_deref().unwrap_or("")))
            .join("\n");
          format!("Length {length}:\n{names}")
        })
        .join("\n\n");

      make_error!("Sequences are expected to all have the same length, but found the following lengths:\n\n{message}")
    },
  }
  .wrap_err("When calculating length of sequences")
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
    },
  }
  .wrap_err("When calculating length of sequences")
}
