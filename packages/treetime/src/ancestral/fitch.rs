use crate::partition::fitch::partition::PartitionFitch;
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::make_internal_report;

pub(crate) fn ancestral_reconstruction_fitch(graph: &Graph, partitions: &mut [PartitionFitch]) -> Result<(), Report> {
  graph.iter_depth_first_preorder_forward(|node| run_fitch_reconstruction(partitions, &node))
}

fn run_fitch_reconstruction(partitions: &mut [PartitionFitch], node: &GraphNodeForward) -> Result<(), Report> {
  for partition in partitions.iter_mut() {
    let alphabet = partition.alphabet.clone();

    let mut sequence = if !node.is_root {
      let (parent, edge) = get_exactly_one(&node.parent_keys)?;
      let mut sequence = partition.nodes[parent].seq.sequence.clone();
      let edge_part = &partition.edges[edge];

      for sub in edge_part.fitch_subs() {
        sequence[sub.pos()] = sub.qry();
      }

      for indel in &edge_part.indels {
        if indel.is_deletion() {
          sequence[indel.range.0..indel.range.1].fill(alphabet.gap());
        } else {
          sequence[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
        }
      }
      sequence
    } else {
      partition.nodes[&node.key].seq.sequence.clone()
    };

    let node_data = partition
      .nodes
      .get_mut(&node.key)
      .ok_or_else(|| make_internal_report!("Node {} has no Fitch data", node.key))?;
    let seq = &mut node_data.seq;

    for r in &mut seq.unknown {
      sequence[r.0..r.1].fill(alphabet.unknown());
    }

    for (pos, states) in &mut seq.fitch.variable {
      sequence[*pos] = alphabet.set_to_char(*states);
    }

    seq.sequence = sequence;
  }
  Ok(())
}
