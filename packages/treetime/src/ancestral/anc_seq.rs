use crate::alphabet::sequence_data::SequenceData;
use crate::ancestral::anc_args::TreetimeAncestralArgs;
use crate::ancestral::anc_graph::{AncestralGraph, NodeType};
use crate::ancestral::run_ancestral_reconstruction::TreetimeAncestralParams;
use crate::io::fasta::FastaRecord;
use itertools::Itertools;

pub fn reconstruct_ancestral_sequences(
  sequence_data: &SequenceData,
  graph: &AncestralGraph,
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
      index: node.key.as_usize(),
      seq_name: payload.name.clone(),
      seq,
    })
  })
}
