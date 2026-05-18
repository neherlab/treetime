use crate::alphabet::alphabet::{FILL_CHAR, NON_CHAR};
use crate::ancestral::fitch_indel::{resolve_indels_backward, resolve_indels_forward};
use crate::ancestral::fitch_sub::{
  finalize_sequence_forward, resolve_fixed_positions_backward, resolve_nonroot_substitutions_forward,
  resolve_root_forward, resolve_variable_positions_backward,
};
use crate::make_report;
use crate::representation::partition::fitch::PartitionFitch;
use crate::representation::partition::traits::PartitionCompressed;
use crate::representation::payload::ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::payload::sparse::{
  FitchSeqDistribution, MarginalSparseSeqDistribution, SparseEdgePartition, SparseNodePartition, SparseSeqInfo,
};
use crate::seq::composition::Composition;
use eyre::{Report, WrapErr};
use itertools::Itertools;
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;
use treetime_graph::breadth_first::GraphTraversalContinuation;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::{GraphNodeBackward, GraphNodeForward};
use treetime_graph::node::{GraphNode, NodeAncestralOps};
use treetime_io::fasta::FastaRecord;
use treetime_primitives::{AlphabetLike, Seq, seq};
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range_intersection::range_intersection_iter;
use treetime_utils::sync::mutex::extract_parallel_error;

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
  for leaf in graph.get_leaves() {
    let leaf_key = leaf.read_arc().key();
    let mut leaf = leaf.read_arc().payload().write_arc();

    let leaf_name = leaf.name().ok_or_else(|| {
      make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
    })?.as_ref().to_owned();

    let leaf_fasta = aln
      .iter()
      .find(|fasta| fasta.seq_name == leaf_name)
      // TODO: we could optionally emit a warning here and continue with a sequence that is missing
      .ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;

    leaf.set_desc(leaf_fasta.desc.clone());

    partitions.iter().try_for_each(|partition| -> Result<(), Report> {
      let mut partition = partition.write_arc();
      let alphabet = &partition.alphabet().clone(); // TODO: avoid clone

      partition
        .nodes_mut()
        .insert(leaf_key, SparseNodePartition::new(&leaf_fasta.seq, alphabet)?);

      Ok(())
    })?;
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

pub(crate) fn fitch_backward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>])
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionCompressed,
{
  graph.par_iter_breadth_first_backward(|node| {
    run_fitch_backward(partitions, &node);
    GraphTraversalContinuation::Continue
  });
}

fn run_fitch_backward<N, E, P>(partitions: &[Arc<RwLock<P>>], node: &GraphNodeBackward<N, E, ()>)
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionCompressed,
{
  if node.is_leaf {
    return;
  }

  for partition in partitions {
    let mut partition = partition.write_arc();

    let children = node
      .child_keys
      .iter()
      .map(|(child, edge)| (&partition.node(child).seq, partition.edge(edge)))
      .collect_vec();

    // determine parts of the sequence that are unknown, gaps in all children
    let mut gaps = range_intersection_iter(children.iter().map(|(c, _)| &c.gaps)).collect_vec();
    let unknown = range_intersection_iter(children.iter().map(|(c, _)| &c.unknown)).collect_vec();
    // non_char are ranges that are either unknown or gaps in all children (can differ from the union of gaps and unknown)
    let non_char = range_intersection_iter(children.iter().map(|(c, _)| &c.non_char)).collect_vec();

    let mut sequence = seq![FILL_CHAR; partition.length()];
    for r in &non_char {
      sequence[r.0..r.1].fill(NON_CHAR);
    }

    // Process all positions where the children are variable.
    // Need to account for parts of the sequence transmitted along edges.
    let mut variable = resolve_variable_positions_backward(&children, &non_char, &mut sequence);
    // Process all positions where the children are fixed or completely unknown in some children.
    resolve_fixed_positions_backward(&children, partition.alphabet(), &mut sequence, &mut variable);

    // Resolve variable indels from child gap disagreements
    let child_gaps_vec: Vec<Vec<(usize, usize)>> = children.iter().map(|(c, _)| c.gaps.clone()).collect_vec();
    let child_variable_indels: Vec<&_> = children.iter().map(|(c, _)| &c.fitch.variable_indel).collect_vec();

    let indels_bw = resolve_indels_backward(&child_gaps_vec, &child_variable_indels, partition.length());

    let new_node_data = SparseNodePartition {
      seq: SparseSeqInfo {
        gaps: {
          gaps.extend(indels_bw.resolved_gaps);
          gaps
        },
        unknown,
        non_char,
        fitch: FitchSeqDistribution {
          variable,
          variable_indel: indels_bw.variable_indel,
          composition: Composition::new(partition.alphabet().chars(), partition.alphabet().gap()),
          chosen_state: btreemap! {},
        },
        sequence,
        composition: Composition::new(partition.alphabet().chars(), partition.alphabet().gap()),
      },
      profile: MarginalSparseSeqDistribution {
        variable: btreemap! {},
        variable_indel: btreemap! {},
        fixed: btreemap! {},
        fixed_counts: Composition::new(partition.alphabet().chars(), partition.alphabet().gap()),
        log_lh: 0.0,
      },
    };

    partition.nodes_mut().insert(node.key, new_node_data);
  }
}

pub(crate) fn fitch_forward<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>]) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionCompressed,
{
  let error: Arc<parking_lot::Mutex<Option<Report>>> = Arc::new(parking_lot::Mutex::new(None));
  graph.par_iter_breadth_first_forward(|node| {
    if let Err(e) = run_fitch_forward(partitions, &node) {
      let mut guard = error.lock();
      if guard.is_none() {
        *guard = Some(e);
      }
      return GraphTraversalContinuation::Stop;
    }
    GraphTraversalContinuation::Continue
  });
  extract_parallel_error(error)
}

fn run_fitch_forward<N, E, P>(partitions: &[Arc<RwLock<P>>], node: &GraphNodeForward<N, E, ()>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionCompressed,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    let alphabet = &partition.alphabet().clone(); // TODO: avoid clone

    let mut node_data = partition.nodes_mut().remove(&node.key).unwrap();

    if node.is_root {
      let seq = &mut node_data.seq;
      resolve_root_forward(
        &mut seq.sequence,
        &mut seq.gaps,
        &seq.fitch.variable,
        &seq.fitch.variable_indel,
        &mut seq.fitch.chosen_state,
      );
    } else {
      let (parent_key, edge_key) =
        get_exactly_one(&node.parent_keys).wrap_err("Multiple parent nodes are not yet supported")?;

      let parent = &partition.node(parent_key).seq;

      let seq = &mut node_data.seq;
      seq.composition = parent.composition.clone();

      // fill in the indeterminate positions by copying the parent (new gaps in the node will be introduced later)
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

      // Resolve variable indels using parent gap context
      let indels = resolve_indels_forward(
        &seq.fitch.variable_indel,
        &mut seq.gaps,
        &parent.gaps,
        &parent.sequence,
        &seq.sequence,
      );
      for indel in &indels {
        seq.composition.add_indel(indel);
      }
      for r in &seq.unknown {
        // this might result in compensating addition/deletions of Ns already present in the parent
        for pos in r.0..r.1 {
          seq.composition.adjust_count(seq.sequence[pos], -1);
        }
        seq
          .composition
          .adjust_count(alphabet.unknown(), r.1 as isize - r.0 as isize);
      }

      let edge = partition.edge_mut(edge_key);
      edge.extend_fitch_subs(subs);
      edge.indels.extend(indels);
    }

    let seq = &mut node_data.seq;
    finalize_sequence_forward(
      &mut seq.sequence,
      &seq.gaps,
      &seq.unknown,
      &mut seq.composition,
      alphabet,
      node.is_root,
    );

    partition.nodes_mut().insert(node.key, node_data);
  }
  Ok(())
}

fn fitch_cleanup<N, E, P>(graph: &Graph<N, E, ()>, partitions: &[Arc<RwLock<P>>])
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionCompressed,
{
  graph.par_iter_breadth_first_forward(|node| run_fitch_forward_cleanup(&node, partitions));
}

fn run_fitch_forward_cleanup<N, E, P>(
  node: &GraphNodeForward<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
) -> GraphTraversalContinuation
where
  N: GraphNode,
  E: GraphEdge,
  P: PartitionCompressed,
{
  for partition in partitions {
    let mut partition = partition.write_arc();
    let seq = &mut partition.node_mut(&node.key).seq;

    // delete the variable position everywhere except of leaves
    if !node.is_leaf {
      seq.fitch.variable = btreemap! {};
    }

    seq.fitch.composition = seq.composition.clone();
    for p in seq.fitch.variable.values() {
      if let Some(state) = p.get_one_maybe() {
        seq.fitch.composition.adjust_count(state, -1);
      }
    }

    // Keep the exact reconstructed sequence on every node. Sparse marginal
    // passes need an authoritative reference state for fixed-site lookups at
    // ambiguous-variable positions.
  }

  GraphTraversalContinuation::Continue
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
  fitch_backward(graph, partitions);
  fitch_forward(graph, partitions)?;
  fitch_cleanup(graph, partitions);
  Ok(())
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
  graph.try_iter_depth_first_preorder_forward(|node| {
    run_fitch_reconstruction(include_leaves, partitions, &mut visitor, &node)
  })
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
