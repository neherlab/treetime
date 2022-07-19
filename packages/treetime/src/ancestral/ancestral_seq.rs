use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::ancestral_graph::{Graph, NodeType};
use crate::ancestral::run_ancestral::TreetimeAncestralParams;
use crate::cli::treetime_cli::TreetimeAncestralArgs;
use crate::io::fasta::FastaRecord;
use itertools::Itertools;

pub fn ancestral_reconstruct_sequences(
  sequence_data: &SequenceData,
  graph: &Graph,
  ancestral_args: &TreetimeAncestralArgs,
  ancestral_params: &TreetimeAncestralParams,
) -> Vec<FastaRecord> {
  graph.filter_map(|node| {
    let payload = &node.payload.read();

    let seq = if let (NodeType::Leaf(name), true) = (&payload.node_type, !ancestral_args.reconstruct_tip_states) {
      sequence_data
        .get_full(name)
        .ok_or_else(|| format!("Sequence not found: '{name}'"))
        .unwrap()
    } else {
      sequence_data.decompress(&payload.seq)
    }
    .iter()
    .join("");

    Some(FastaRecord {
      index: node.key,
      seq_name: payload.name.clone(),
      seq,
    })
  })
}
