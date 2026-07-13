use crate::alphabet::alphabet::{Alphabet, FILL_CHAR, NON_CHAR};
use crate::ancestral::fitch_indel::{compute_node_ranges, resolve_indels_backward, resolve_indels_forward};
use crate::ancestral::fitch_sub::{
  finalize_sequence_forward, resolve_fixed_positions_backward, resolve_nonroot_substitutions_forward,
  resolve_root_forward, resolve_variable_positions_backward,
};
use crate::make_report;
use crate::partition::fitch::PartitionFitch;
use crate::partition::indexed_pass::{IndexedPass, IndexedPassSlot};
use crate::partition::sparse::{
  FitchSeqDistribution, SparseEdgePartition, SparseNodePartition, SparseSeqDistribution, SparseSeqInfo,
};
use crate::partition::traits::PartitionCompressed;
use crate::payload::ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::seq::alignment::get_common_length;
use crate::seq::composition::Composition;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use parking_lot::RwLock;
use rayon::prelude::*;
use std::collections::{BTreeMap, BTreeSet};
use std::sync::Arc;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::{GraphNode, NodeAncestralOps};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{AlphabetLike, Seq, seq};
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::sync::mutex::unwrap_arc_rwlock;

pub fn create_fitch_partition<N, E>(
  graph: &Graph<N, E, ()>,
  index: usize,
  alphabet: Alphabet,
  aln: &[FastaRecord],
) -> Result<PartitionFitch, Report>
where
  N: NodeAncestralOps,
  E: GraphEdge,
{
  let length = get_common_length(aln)?;
  let partition = Arc::new(RwLock::new(PartitionFitch {
    index,
    alphabet,
    length,
    nodes: btreemap! {},
    edges: btreemap! {},
  }));
  compress_sequences(graph, std::slice::from_ref(&partition), aln)?;
  unwrap_arc_rwlock(partition)
}

pub(crate) fn attach_seqs_to_graph<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  aln: &[FastaRecord],
) -> Result<(), Report>
where
  N: NodeAncestralOps,
  E: GraphEdge,
  P: PartitionCompressed,
{
  let aln_by_name = aln.iter().fold(BTreeMap::new(), |mut records, record| {
    records.entry(record.seq_name.as_str()).or_insert(record);
    records
  });
  let leaf_records = graph
    .get_leaves()
    .into_par_iter()
    .map(|leaf| -> Result<_, Report> {
      let leaf = leaf.read_arc();
      let leaf_key = leaf.key();
      let mut leaf_payload = leaf.payload().write_arc();
      let leaf_name = leaf_payload
        .name()
        .ok_or_else(|| {
          make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
        })?
        .as_ref()
        .to_owned();
      let leaf_fasta = aln_by_name
        .get(leaf_name.as_str())
        .copied()
        // Every leaf has a sequence record after alignment completion.
        .ok_or_else(|| make_report!("Leaf sequence not found after alignment completion: '{leaf_name}'"))?;
      leaf_payload.set_desc(leaf_fasta.desc.clone());
      Ok((leaf_key, leaf_fasta))
    })
    .collect::<Result<Vec<_>, Report>>()?;

  for partition in partitions {
    let alphabet = partition.read_arc().alphabet().clone();
    let nodes = leaf_records
      .par_iter()
      .map(|(leaf_key, leaf_fasta)| SparseNodePartition::new(&leaf_fasta.seq, &alphabet).map(|node| (*leaf_key, node)))
      .collect::<Result<BTreeMap<_, _>, Report>>()?;
    partition.write_arc().nodes_mut().extend(nodes);
  }

  for edge in graph.get_edges() {
    let edge_key = edge.read_arc().key();
    partitions.iter().try_for_each(|partition| -> Result<(), Report> {
      let mut partition = partition.write_arc();
      partition.edges_mut().insert(edge_key, SparseEdgePartition::default());
      Ok(())
    })?;
  }

  Ok(())
}

pub(crate) fn fitch_backward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionCompressed,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    let alphabet = partition.alphabet().clone();
    let length = partition.length();
    let (nodes, edges) = partition.storage_mut();
    let mut pass = IndexedPass::new(graph, nodes, edges, |_| Ok(SparseNodePartition::empty(&alphabet)))?;
    let result = pass.try_for_each_backward_frontier(|node_indices, _, _, completed, frontier| {
      frontier
        .par_iter_mut()
        .try_for_each(|slot| run_fitch_backward_indexed(graph, &alphabet, length, node_indices, completed, slot))
    });
    let (nodes, edges) = pass.into_maps()?;
    *partition.nodes_mut() = nodes;
    *partition.edges_mut() = edges;
    result?;
  }
  Ok(())
}

fn run_fitch_backward_indexed<N, E>(
  graph: &Graph<N, E, ()>,
  alphabet: &Alphabet,
  length: usize,
  node_indices: &[Option<usize>],
  completed: &[IndexedPassSlot<SparseNodePartition, SparseEdgePartition>],
  slot: &mut IndexedPassSlot<SparseNodePartition, SparseEdgePartition>,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let graph_node = graph.get_node(slot.key).expect("Indexed node must exist in graph");
  let graph_node = graph_node.read_arc();
  if graph_node.is_leaf() {
    return Ok(());
  }
  let child_keys = graph.children_of(&graph_node);
  let children = child_keys
    .iter()
    .map(|(child, _)| {
      let child_key = child.read_arc().key();
      let child_index = node_indices[child_key.as_usize()].expect("Indexed child must have a slot");
      let child = &completed[child_index];
      let edge = &child
        .parent_edge
        .as_ref()
        .expect("Non-root indexed node must own its parent edge")
        .1;
      (&child.node.seq, edge)
    })
    .collect_vec();

  let child_non_chars: Vec<&Vec<(usize, usize)>> = children.iter().map(|(c, _)| &c.non_char).collect_vec();
  let child_gaps: Vec<&Vec<(usize, usize)>> = children.iter().map(|(c, _)| &c.gaps).collect_vec();

  let ranges = compute_node_ranges(&child_non_chars, &child_gaps);
  let non_char = ranges.non_char;
  let unknown = ranges.unknown;

  let mut sequence = seq![FILL_CHAR; length];
  for r in &non_char {
    sequence[r.0..r.1].fill(NON_CHAR);
  }

  let mut variable = resolve_variable_positions_backward(&children, &non_char, &mut sequence);
  resolve_fixed_positions_backward(&children, alphabet, &mut sequence, &mut variable);

  let child_unknown: Vec<&Vec<(usize, usize)>> = children.iter().map(|(c, _)| &c.unknown).collect_vec();
  let child_variable_indels: Vec<&_> = children.iter().map(|(c, _)| &c.fitch.variable_indel).collect_vec();

  let indels_bw = resolve_indels_backward(&child_gaps, &child_unknown, &child_variable_indels, length);

  slot.node = SparseNodePartition {
    seq: SparseSeqInfo {
      gaps: indels_bw.resolved_gaps,
      unknown,
      non_char,
      fitch: FitchSeqDistribution {
        variable,
        variable_indel: indels_bw.variable_indel,
        chosen_state: btreemap! {},
      },
      sequence,
      composition: Composition::new(alphabet.chars(), alphabet.gap()),
    },
    profile: SparseSeqDistribution {
      variable: btreemap! {},
      variable_indel: BTreeSet::new(),
      fixed: btreemap! {},
      fixed_counts: Composition::new(alphabet.chars(), alphabet.gap()),
      log_lh: 0.0,
    },
  };
  Ok(())
}

pub(crate) fn fitch_forward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionCompressed,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    let alphabet = partition.alphabet().clone();
    let (nodes, edges) = partition.storage_mut();
    let mut pass = IndexedPass::new(graph, nodes, edges, |key| {
      Err(make_report!(
        "Partition node {key} is missing before the Fitch forward pass"
      ))
    })?;
    let result = pass.try_for_each_forward_frontier(|node_indices, _, _, completed_start, frontier, completed| {
      frontier
        .par_iter_mut()
        .try_for_each(|slot| run_fitch_forward_indexed(&alphabet, node_indices, completed_start, completed, slot))
    });
    let (nodes, edges) = pass.into_maps()?;
    *partition.nodes_mut() = nodes;
    *partition.edges_mut() = edges;
    result?;
  }
  Ok(())
}

fn run_fitch_forward_indexed(
  alphabet: &Alphabet,
  node_indices: &[Option<usize>],
  completed_start: usize,
  completed: &[IndexedPassSlot<SparseNodePartition, SparseEdgePartition>],
  slot: &mut IndexedPassSlot<SparseNodePartition, SparseEdgePartition>,
) -> Result<(), Report> {
  let node_data = &mut slot.node;
  if let Some(parent_key) = slot.parent_key {
    let (_, edge) = slot
      .parent_edge
      .as_mut()
      .expect("Non-root node must own its parent edge");
    let parent_index = node_indices[parent_key.as_usize()].expect("Indexed parent must have a slot");
    let parent = &completed[parent_index - completed_start].node.seq;
    let seq = &mut node_data.seq;
    seq.composition = parent.composition.clone();

    for r in &seq.non_char {
      seq.sequence[r.0..r.1].clone_from_slice(&parent.sequence[r.0..r.1]);
    }

    let subs = resolve_nonroot_substitutions_forward(
      &mut seq.sequence,
      &seq.gaps,
      &mut seq.fitch.variable,
      &mut seq.fitch.chosen_state,
      &mut seq.composition,
      parent,
      alphabet,
    )?;

    let (indels, new_gaps) = resolve_indels_forward(
      &seq.fitch.variable_indel,
      &seq.gaps,
      &seq.non_char,
      &parent.gaps,
      &parent.sequence,
      &seq.sequence,
    );
    seq.gaps = new_gaps;
    for indel in &indels {
      seq.composition.add_indel(indel);
    }
    for r in &seq.unknown {
      for pos in r.0..r.1 {
        seq.composition.adjust_count(seq.sequence[pos], -1);
      }
      seq
        .composition
        .adjust_count(alphabet.unknown(), r.1 as isize - r.0 as isize);
    }

    edge.extend_fitch_subs(subs);
    edge.indels.extend(indels);
  } else {
    let seq = &mut node_data.seq;
    resolve_root_forward(&mut seq.sequence, &seq.fitch.variable, &mut seq.fitch.chosen_state);
  }

  let seq = &mut node_data.seq;
  finalize_sequence_forward(
    &mut seq.sequence,
    &seq.gaps,
    &seq.unknown,
    &mut seq.composition,
    alphabet,
    slot.parent_key.is_none(),
  );
  Ok(())
}

fn fitch_cleanup<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionCompressed,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    for (key, node) in partition.nodes_mut() {
      if !graph.is_leaf(*key) {
        node.seq.fitch.variable = btreemap! {};
      }
    }
  }
  Ok(())
}

pub fn compress_sequences<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
  aln: &[FastaRecord],
) -> Result<(), Report>
where
  N: NodeAncestralOps,
  E: GraphEdge,
  P: PartitionCompressed,
{
  attach_seqs_to_graph(graph, partitions, aln)?;
  fitch_backward(graph, partitions)?;
  fitch_forward(graph, partitions)?;
  fitch_cleanup(graph, partitions)
}

/// Reconstruct ancestral sequences using Fitch parsimony.
///
/// Calls visitor function for every ancestral node, providing the node itself and its reconstructed sequence.
/// Optionally reconstructs leaf sequences.
pub fn ancestral_reconstruction_fitch(
  graph: &GraphAncestral,
  include_leaves: bool,
  partitions: &[Arc<RwLock<PartitionFitch>>],
  mut visitor: impl FnMut(&GraphNodeForward<NodeAncestral, EdgeAncestral, ()>, &Seq) -> Result<(), Report>,
) -> Result<(), Report> {
  graph
    .iter_depth_first_preorder_forward(|node| run_fitch_reconstruction(include_leaves, partitions, &mut visitor, &node))
}

fn run_fitch_reconstruction(
  include_leaves: bool,
  partitions: &[Arc<RwLock<PartitionFitch>>],
  mut visitor: impl FnMut(&GraphNodeForward<NodeAncestral, EdgeAncestral, ()>, &Seq) -> Result<(), Report>,
  node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  if !include_leaves && node.is_leaf {
    return Ok(());
  }

  for partition in partitions {
    let alphabet = &partition.read_arc().alphabet.clone(); // TODO: avoid clone

    let mut sequence = if !node.is_root {
      let partition = partition.read_arc();
      let (parent, edge) = get_exactly_one(&node.parent_keys).unwrap();
      let mut sequence = partition.nodes[parent].seq.sequence.clone();
      let edge_part = &partition.edges[edge];

      for sub in edge_part.fitch_subs() {
        sequence[sub.pos()] = sub.qry();
      }

      for indel in &edge_part.indels {
        if indel.deletion {
          sequence[indel.range.0..indel.range.1].fill(alphabet.gap());
        } else {
          sequence[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
        }
      }
      sequence
    } else {
      let partition = partition.read_arc();
      partition.nodes[&node.key].seq.sequence.clone()
    };

    let mut partition = partition.write_arc();
    let node_data = partition.nodes.get_mut(&node.key).unwrap();
    let seq = &mut node_data.seq;

    for r in &mut seq.unknown {
      sequence[r.0..r.1].fill(alphabet.unknown());
    }

    for (pos, states) in &mut seq.fitch.variable {
      sequence[*pos] = alphabet.set_to_char(*states);
    }

    seq.sequence = sequence;

    visitor(node, &seq.sequence)?;
  }
  Ok(())
}
