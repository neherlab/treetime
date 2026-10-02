use crate::partition::fitch::partition::PartitionFitch;
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_utils::collections::container::get_exactly_one;

pub(crate) fn ancestral_reconstruction_fitch(
  graph: &Graph,
  include_leaves: bool,
  partitions: &mut [PartitionFitch],
) -> Result<Vec<GraphNodeKey>, Report> {
  let mut emitted_nodes = Vec::new();
  graph.iter_depth_first_preorder_forward(|node| {
    if run_fitch_reconstruction(include_leaves, partitions, &node)? {
      emitted_nodes.push(node.key);
    }
    Ok(())
  })?;
  Ok(emitted_nodes)
}

#[allow(
  clippy::unwrap_used,
  reason = "unwrap on a value an upstream invariant guarantees is present"
)]
fn run_fitch_reconstruction(
  include_leaves: bool,
  partitions: &mut [PartitionFitch],
  node: &GraphNodeForward,
) -> Result<bool, Report> {
  if !include_leaves && node.is_leaf {
    return Ok(false);
  }

  for partition in partitions.iter_mut() {
    let alphabet = partition.alphabet.clone();

    let mut sequence = if !node.is_root {
      let (parent, edge) = get_exactly_one(&node.parent_keys).unwrap();
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

    let node_data = partition.nodes.get_mut(&node.key).unwrap();
    let seq = &mut node_data.seq;

    for r in &mut seq.unknown {
      sequence[r.0..r.1].fill(alphabet.unknown());
    }

    for (pos, states) in &mut seq.fitch.variable {
      sequence[*pos] = alphabet.set_to_char(*states);
    }

    seq.sequence = sequence;
  }
  Ok(true)
}
