use crate::alphabet::alphabet::{Alphabet, FILL_CHAR, NON_CHAR};
use crate::ancestral::fitch_indel::{compute_node_ranges, resolve_indels_backward, resolve_indels_forward};
use crate::ancestral::fitch_sub::{
  discover_fixed_disagreements_backward, finalize_sequence_forward, resolve_nonroot_substitutions_forward,
  resolve_root_forward, resolve_variable_positions_backward,
};
use crate::make_report;
use crate::partition::fitch::partition::PartitionFitch;
use crate::partition::storage::sparse::{FitchNodeData, FitchSeqDistribution, FitchSeqInfo, SparseEdgeObs};
use crate::seq::alignment::NodeSeqInput;
use crate::seq::alignment::get_common_length_of_node_inputs;
use crate::seq::composition::Composition;
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use rayon::prelude::*;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::graph_traverse::GraphNodeForward;
use treetime_graph::node::GraphNodeKey;
use treetime_graph::pass::{GraphPass, GraphPassBackwardContext, GraphPassForwardContext, GraphPassNodeOutput};
use treetime_primitives::{AlphabetLike, seq};
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range_union::range_union;

pub fn create_fitch_partition(
  graph: &Graph,
  index: usize,
  alphabet: Alphabet,
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
) -> Result<PartitionFitch, Report> {
  let length = get_common_length_of_node_inputs(node_inputs)?;
  let mut partition = PartitionFitch {
    index,
    alphabet,
    length,
    nodes: btreemap! {},
    edges: btreemap! {},
  };
  compress_sequences(graph, &mut partition, node_inputs)?;
  Ok(partition)
}

pub(crate) fn attach_seqs_to_graph(
  graph: &Graph,
  partition: &mut PartitionFitch,
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
) -> Result<(), Report> {
  let leaf_records = graph
    .get_leaves()
    .collect::<Vec<_>>()
    .into_par_iter()
    .map(|leaf| -> Result<_, Report> {
      let leaf_key = leaf.key();
      let node = &node_inputs[&leaf_key];
      let seq = node
        .seq
        .as_ref()
        // Every leaf has a sequence after alignment completion.
        .ok_or_else(|| {
          make_report!(
            "Leaf sequence not found after alignment completion: '{}'",
            node.name.as_deref().unwrap_or("")
          )
        })?;
      Ok((leaf_key, seq))
    })
    .collect::<Result<Vec<_>, Report>>()?;

  let alphabet = partition.alphabet.clone();
  let nodes = leaf_records
    .par_iter()
    .map(|(leaf_key, seq)| FitchNodeData::new(seq, &alphabet).map(|node| (*leaf_key, node)))
    .collect::<Result<BTreeMap<_, _>, Report>>()?;
  partition.nodes.extend(nodes);

  for edge in graph.get_edges() {
    let edge_key = edge.key();
    partition.edges.insert(edge_key, SparseEdgeObs::default());
  }

  Ok(())
}

pub(crate) fn fitch_backward(graph: &Graph, partition: &mut PartitionFitch) -> Result<(), Report> {
  let alphabet = partition.alphabet.clone();
  let length = partition.length;
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_backward(
    &partition.nodes,
    &partition.edges,
    |_| Ok(FitchNodeData::empty(&alphabet)),
    |context| run_fitch_backward_indexed(&alphabet, length, &context),
  )?;
  partition.nodes = outputs.nodes;
  partition.edges = outputs.edges;
  Ok(())
}

#[allow(clippy::expect_used, reason = "expect on a value an upstream invariant guarantees is present")]
fn run_fitch_backward_indexed(
  alphabet: &Alphabet,
  length: usize,
  context: &GraphPassBackwardContext<'_, FitchNodeData, SparseEdgeObs, FitchNodeData, SparseEdgeObs>,
) -> Result<GraphPassNodeOutput<FitchNodeData, SparseEdgeObs>, Report> {
  if context.is_leaf {
    // A leaf keeps its attached Fitch data unchanged and returns its parent edge untouched so the
    // forward pass keeps the edge entries it depends on. A single-node tree, where the leaf is also the
    // root, has no parent edge and returns no message.
    let node = context.input.clone();
    let parent_message = context.parent_edge.map(|(_, edge)| edge.clone());
    return Ok(GraphPassNodeOutput { node, parent_message });
  }

  // Children arrive in the graph's canonical `children_of` order, so every child is fetched and folded
  // in that order, keeping the parsimony result byte-for-byte identical.
  let children = context
    .children
    .iter()
    .map(|child| {
      let edge_data = child
        .edge
        .expect("Backward child edge message must be published before its parent");
      (&child.node.seq, edge_data)
    })
    .collect_vec();

  let child_non_chars: Vec<&Vec<(usize, usize)>> = children.iter().map(|(c, _)| &c.non_char).collect_vec();
  let child_gaps: Vec<&Vec<(usize, usize)>> = children.iter().map(|(c, _)| &c.gaps).collect_vec();

  let ranges = compute_node_ranges(&child_non_chars, &child_gaps);
  let unknown = ranges.unknown;

  let child_unknown: Vec<&Vec<(usize, usize)>> = children.iter().map(|(c, _)| &c.unknown).collect_vec();
  let child_variable_indels: Vec<&_> = children.iter().map(|(c, _)| &c.fitch.variable_indel).collect_vec();

  let indels_bw = resolve_indels_backward(&child_gaps, &child_unknown, &child_variable_indels, length);

  // A resolved gap is a position with no character state, so it has to be masked like every other
  // non-char position. `compute_node_ranges` intersects the children's `non_char`, which drops any
  // column a child left as `variable_indel`, while `resolve_indels_backward` still resolves such a
  // column to a gap (it counts `variable_indel` as gap-compatible). Taking the union keeps `gaps` a
  // subset of `non_char`, the invariant leaves (`FitchNodeData::new`) and the dense
  // representation (`DenseSeqInfo::new`) already hold. Without it a single determined residue
  // stranded inside a missing-data run keeps a character state at a position the node reports as
  // deleted, and the forward pass then emits a substitution inside its own deletion.
  let non_char = range_union(&[ranges.non_char, indels_bw.resolved_gaps.clone()]);

  let mut sequence = seq![FILL_CHAR; length];
  for r in &non_char {
    sequence[r.0..r.1].fill(NON_CHAR);
  }

  // Discovery first, resolution second. The discovery pass only flags positions where children
  // hold differing canonical states; the resolution pass then recomputes every candidate position
  // from all children, so each child is counted exactly once. Running them the other way round
  // let both passes fold child states into the same map, which is harmless for a union but not
  // for the plurality rule.
  let discovered = discover_fixed_disagreements_backward(&children, alphabet, &mut sequence);
  let variable = resolve_variable_positions_backward(&children, &discovered, &non_char, &mut sequence);

  let node = FitchNodeData {
    seq: FitchSeqInfo {
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
  };

  // Fitch backward computes only node data. A non-root node returns its parent edge unchanged so the
  // forward pass keeps its edge entries; the root has no parent edge.
  let parent_message = context.parent_edge.map(|(_, edge)| edge.clone());
  Ok(GraphPassNodeOutput { node, parent_message })
}

pub(crate) fn fitch_forward(graph: &Graph, partition: &mut PartitionFitch) -> Result<(), Report> {
  let alphabet = partition.alphabet.clone();
  let pass = GraphPass::new(graph)?;
  let outputs = pass.map_forward(
    &partition.nodes,
    &partition.edges,
    |key| {
      Err(make_report!(
        "Partition node {key} is missing before the Fitch forward pass"
      ))
    },
    |context| run_fitch_forward_indexed(&alphabet, &context),
  )?;
  partition.nodes = outputs.nodes;
  partition.edges = outputs.edges;
  Ok(())
}

#[allow(clippy::as_conversions, clippy::expect_used, reason = "count/index numeric cast is exact for the domain range; expect on a value an upstream invariant guarantees is present")]
fn run_fitch_forward_indexed(
  alphabet: &Alphabet,
  context: &GraphPassForwardContext<'_, FitchNodeData, SparseEdgeObs, FitchNodeData>,
) -> Result<GraphPassNodeOutput<FitchNodeData, SparseEdgeObs>, Report> {
  let mut node = context.input.clone();

  // The forward pass produces the durable edge data for the asymmetric Fitch case, so a non-root
  // node clones and extends its parent edge rather than building a fresh one: the edge arrived through
  // the backward pass carrying fields that must survive here.
  let parent_message = if let Some((_, edge)) = context.parent_edge {
    let mut edge = edge.clone();
    let parent = &context.parent.expect("Non-root node must have a parent").seq;
    let seq = &mut node.seq;
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
    // The forward pass widens `gaps` with gaps inherited from the parent, so re-establish
    // `gaps ⊆ non_char` here too. Must run after `resolve_indels_forward`, which distinguishes
    // insertions by testing `non_char` as it stood during the backward pass.
    seq.non_char = range_union(&[seq.non_char.clone(), seq.gaps.clone()]);
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
    Some(edge)
  } else {
    let seq = &mut node.seq;
    resolve_root_forward(
      &mut seq.sequence,
      &seq.fitch.variable,
      &mut seq.fitch.chosen_state,
      alphabet,
    );
    None
  };

  let seq = &mut node.seq;
  finalize_sequence_forward(
    &mut seq.sequence,
    &seq.gaps,
    &seq.unknown,
    &mut seq.composition,
    alphabet,
    context.is_root,
  );

  Ok(GraphPassNodeOutput { node, parent_message })
}

fn fitch_cleanup(graph: &Graph, partition: &mut PartitionFitch) -> Result<(), Report> {
  for (key, node) in &mut partition.nodes {
    if !graph.is_leaf(*key) {
      node.seq.fitch.variable = btreemap! {};
    }
  }
  Ok(())
}

pub fn compress_sequences(
  graph: &Graph,
  partition: &mut PartitionFitch,
  node_inputs: &BTreeMap<GraphNodeKey, NodeSeqInput>,
) -> Result<(), Report> {
  attach_seqs_to_graph(graph, partition, node_inputs)?;
  fitch_backward(graph, partition)?;
  fitch_forward(graph, partition)?;
  fitch_cleanup(graph, partition)
}

/// Reconstruct ancestral sequences using Fitch parsimony.
///
/// Calls the visitor for every reconstructed node, providing the node itself and its reconstructed
/// sequence, and returns the reconstructed sequences keyed by node id. The returned map is the
/// reconstruction result as a value the command captures; each sequence also stays written into the
/// partition (read by the node-data serializer until the tree writers read the map directly).
/// Optionally reconstructs leaf sequences.
/// Walk the tree in depth-first preorder, writing each node's reconstructed sequence into its partition
/// state (the parent-to-child propagation store the child reads), and return the node ids that emit a
/// sequence, in walk order.
pub fn ancestral_reconstruction_fitch(
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

#[allow(clippy::unwrap_used, reason = "unwrap on a value an upstream invariant guarantees is present")]
/// Reconstruct one node's sequence into its partition state. Returns `true` when the node emits a
/// sequence, `false` for a suppressed tip.
fn run_fitch_reconstruction(
  include_leaves: bool,
  partitions: &mut [PartitionFitch],
  node: &GraphNodeForward,
) -> Result<bool, Report> {
  if !include_leaves && node.is_leaf {
    return Ok(false);
  }

  for partition in partitions.iter_mut() {
    let alphabet = partition.alphabet.clone(); // TODO: avoid clone

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
