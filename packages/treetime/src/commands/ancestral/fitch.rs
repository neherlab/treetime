use crate::alphabet::alphabet::{FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
use crate::representation::partition::fitch::PartitionFitch;
use crate::representation::partition::traits::PartitionCompressed;
use crate::representation::payload::ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::payload::sparse::{
  Deletion, FitchSeqDistribution, MarginalSparseSeqDistribution, SparseEdgePartition, SparseNodePartition,
  SparseSeqInfo,
};
use crate::seq::composition::Composition;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use crate::{make_error, make_report};
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
use treetime_primitives::AlphabetLike;
use treetime_primitives::{BitSet128, Seq, StateSet, StateSetStatus, seq, stateset};
use treetime_utils::collections::container::get_exactly_one;
use treetime_utils::interval::range::range_contains;
use treetime_utils::interval::range_complement::range_complement;
use treetime_utils::interval::range_difference::range_difference;
use treetime_utils::interval::range_intersection::{range_intersection, range_intersection_iter};

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

    let n_children = node.child_keys.len();

    let children = node
      .child_keys
      .iter()
      .map(|(child, edge)| (&partition.node(child).seq, partition.edge(edge)))
      .collect_vec();

    // determine parts of the sequence that are unknown, gaps in all children
    let mut gaps = range_intersection_iter(children.iter().map(|(c, _)| &c.gaps)).collect_vec();
    let unknown = range_intersection_iter(children.iter().map(|(c, _)| &c.unknown)).collect_vec();
    // non_char are ranges that are either unknown or gaps in all children (note that can be different from the union of gaps and unknown)
    let non_char = range_intersection_iter(children.iter().map(|(c, _)| &c.non_char)).collect_vec();
    // calculate the complement of gaps for later look-up
    let non_gap = range_complement(&[(0, partition.length())], &[gaps.clone()]); // FIXME(perf): unnecessary clone

    // what follows could be a function that returns `sequence` and `variable`, takes as arguments children, non_char, alphabet
    let mut sequence = seq![FILL_CHAR; partition.length()];
    for r in &non_char {
      sequence[r.0..r.1].fill(NON_CHAR);
    }

    // Process all positions where the children are variable.
    // Need to account for parts of the sequence transmitted along edges.
    let variable_positions = children
      .iter()
      .flat_map(|(c, _)| c.fitch.variable.keys().copied())
      .unique()
      .collect_vec();

    // Initialization of target data structure (could be done later)
    let mut seq_dis = FitchSeqDistribution {
      variable: btreemap! {},
      variable_indel: btreemap! {},
      composition: Composition::new(partition.alphabet().chars(), partition.alphabet().gap()),
    };

    for pos in variable_positions {
      // Collect child profiles (1D vectors)
      let child_profiles = children
        .iter()
        .filter_map(|(child, edge)| {
          if let Some(transmission) = &edge.transmission {
            if range_contains(transmission, pos) {
              return None; // transmission field is not currently used
            }
          }
          if range_contains(&child.non_char, pos) {
            return None; // this position does not have character state information
          }
          let state = match child.fitch.variable.get(&pos) {
            Some(var_pos) => *var_pos,
            None => StateSet::from_char(child.sequence[pos]),
          };
          Some(state)
        })
        .collect_vec();

      // Calculate Fitch parsimony.
      // If we save the states of the children for each position that is variable in the node,
      // then we would not need the full sequences in the forward pass.
      let intersection = StateSet::from_intersection(&child_profiles);

      match intersection.get() {
        StateSetStatus::Unambiguous(state) => {
          // intersection has a single state, write it
          sequence[pos] = state;
        },
        StateSetStatus::Ambiguous(_) => {
          // more than one possible states
          seq_dis.variable.insert(pos, intersection);
          sequence[pos] = VARIABLE_CHAR;
        },
        StateSetStatus::Empty => {
          let union = StateSet::from_union(&child_profiles);
          seq_dis.variable.insert(pos, union);
          sequence[pos] = VARIABLE_CHAR;
        },
      }
    }

    // Process all positions where the children are fixed or completely unknown in some children.
    for &(child, _) in &children {
      for (pos, parent_state) in sequence.iter_mut().enumerate() {
        let child_state = child.sequence[pos];
        if *parent_state == child_state || *parent_state == NON_CHAR {
          continue; // if parent is equal to child state or we know it's a non-char, skip
        }
        if partition.alphabet().is_canonical(child_state) {
          if *parent_state == FILL_CHAR {
            // if child state is canonical and parent is still FILL_CHAR, set parent_state
            *parent_state = child_state;
          } else {
            // otherwise set or update the variable state
            *seq_dis.variable.entry(pos).or_insert_with(|| stateset! {*parent_state}) += child_state;
            *parent_state = VARIABLE_CHAR;
          }
        }
      }
    }

    // Process insertions and deletions. This also could be a function that returns variable_indel

    // 1) seq_info.gaps is the intersection of child gaps, i.e. this is gap if and only if all children have a gap
    //    --> hence we find positions where children differ in terms of gap presence absence by intersecting
    //        the child gaps with the complement of the parent gaps
    for (child, _) in &children {
      for r in range_intersection(&[non_gap.clone(), child.gaps.clone()]) {
        let indel = seq_dis.variable_indel.entry(r).or_insert_with(|| Deletion {
          deleted: 0,
          present: n_children,
        });
        indel.deleted += 1;
        indel.present -= 1;
      }
    }

    // 2) if a gap is variable in a child and the parent, we need to add this child to the deletion count of this range
    for (child, _) in &children {
      for r in child.fitch.variable_indel.keys() {
        if let Some(indel) = seq_dis.variable_indel.get_mut(r) {
          indel.deleted += 1;
          indel.present -= 1;
        }
      }
    }

    // 3) if all children are compatible with a gap, we add the gap back to the gap collection and remove the
    // variable site (nothing needs doing in the case where all children are compatible with non-gap)
    // this could for example happen if a gap position is variable in a child and thus not part of child.gaps
    seq_dis.variable_indel.retain(|r, indel| {
      if indel.deleted == n_children {
        gaps.push(*r);
        false
      } else {
        true
      }
    });

    let new_node_data = SparseNodePartition {
      seq: SparseSeqInfo {
        gaps,
        unknown,
        non_char,
        fitch: seq_dis,
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

pub(crate) fn fitch_forward<N, E, P>(
  graph: &Graph<N, E, ()>,
  partitions: &[Arc<RwLock<P>>],
) -> Result<(), Report>
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
  match Arc::try_unwrap(error) {
    Ok(mutex) => {
      if let Some(e) = mutex.into_inner() {
        return Err(e);
      }
    }
    Err(arc) => {
      if let Some(e) = arc.lock().take() {
        return Err(e);
      }
    }
  }
  Ok(())
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
      let SparseSeqInfo {
        gaps,
        sequence,
        fitch: FitchSeqDistribution {
          variable,
          variable_indel,
          ..
        },
        ..
      } = &mut node_data.seq;

      for (pos, states) in variable {
        sequence[*pos] = states.get_one();
      }
      // process indels as majority rule at the root
      for (r, indel) in variable_indel.iter() {
        if indel.deleted > indel.present {
          gaps.push(*r);
        }
      }
    } else {
      let SparseSeqInfo {
        gaps,
        unknown,
        sequence,
        composition,
        non_char,
        fitch: FitchSeqDistribution {
          variable,
          variable_indel,
          ..
        },
      } = &mut node_data.seq;

      let (parent_key, edge_key) =
        get_exactly_one(&node.parent_keys).wrap_err("Multiple parent nodes are not yet supported")?;

      let mut subs = vec![];
      let mut indels = vec![];

      let parent = &partition.node(parent_key).seq;
      *composition = parent.composition.clone();

      // fill in the indeterminate positions by copying the parent (note that new gaps in the node will be introduced later)
      for r in non_char {
        sequence[r.0..r.1].clone_from_slice(&parent.sequence[r.0..r.1]);
      }

      // the following two loops modify the sequence, composition, and edge in place and process variable position.
      // for each variable position, pick a state or a mutation
      for (pos, states) in variable.iter_mut() {
        let pnuc = parent.sequence[*pos];
        if alphabet.is_canonical(pnuc) {
          // check whether parent is in child profile (sum>0 --> parent state is in profile)
          if states.contains(pnuc) {
            sequence[*pos] = pnuc;
          } else {
            let cnuc = states.get_one();
            sequence[*pos] = cnuc;
            let m = Sub::new(pnuc, *pos, cnuc)?;
            m.check_determined(alphabet)?;
            composition.add_sub(&m);
            subs.push(m);
          }
        } else if alphabet.is_gap(pnuc) && !range_contains(gaps, *pos) {
          // if parent is gap, but child isn't, we need to resolve variable states
          sequence[*pos] = states.get_one();
        }
      }

      for &pos in parent.fitch.variable.keys() {
        if variable.contains_key(&pos) || range_contains(&parent.gaps, pos) {
          continue;
        }

        // NOTE: access to full_seq would not be necessary if we had saved the
        // child state of variable positions in the backward pass
        let node_nuc = sequence[pos];
        if alphabet.is_canonical(node_nuc) && parent.sequence[pos] != node_nuc {
          let m = Sub::new(parent.sequence[pos], pos, node_nuc)?;
          m.check_determined(alphabet)?;
          composition.add_sub(&m);
          subs.push(m);
        }
      }

      // The following deals with indels.
      // Process indels. Gaps where the children disagree, need to be decided by also looking at parent.
      for (r, indel) in variable_indel.iter() {
        let gap_in_parent = if parent.gaps.contains(r) { 1 } else { 0 };
        if indel.deleted + gap_in_parent > indel.present {
          gaps.push(*r);
          if gap_in_parent == 0 {
            // If the gap is not in parent, add deletion.
            // the sequence that is deleted is the sequence of the parent
            let indel = InDel::del(*r, &parent.sequence[r.0..r.1]);
            composition.add_indel(&indel);
            indels.push(indel);
          }
        } else if gap_in_parent > 0 {
          // Add insertion if gap is present in parent.
          let indel = InDel::ins(*r, &sequence[r.0..r.1]);
          composition.add_indel(&indel);
          indels.push(indel);
        }
      }

      // Process consensus gaps in the node that are not in the parent (deletions)
      for r in range_difference(gaps, &parent.gaps) {
        if variable_indel.contains_key(&r) {
          // all gaps in variable_indel are already processed
          continue;
        }
        let indel = InDel::del(r, &sequence[r.0..r.1]);
        composition.add_indel(&indel);
        indels.push(indel);
      }

      // Process gaps in the parent that are not in the node (insertions)
      for r in range_difference(&parent.gaps, gaps) {
        if variable_indel.contains_key(&r) {
          // all gaps in variable_indel are already processed
          continue;
        }
        let indel = InDel::ins(r, &sequence[r.0..r.1]);
        composition.add_indel(&indel);
        indels.push(indel);
      }
      for r in unknown.iter() {
        // this might result in compensating addition/deletions of Ns already present in the parent
        for pos in r.0..r.1 {
          composition.adjust_count(sequence[pos], -1);
        }
        composition.adjust_count(alphabet.unknown(), r.1 as isize - r.0 as isize);
      }

      {
        let edge = partition.edge_mut(edge_key);
        edge.subs.extend(subs);
        edge.indels.extend(indels);
      }
    }

    let SparseSeqInfo {
      gaps,
      unknown,
      sequence,
      composition,
      ..
    } = &mut node_data.seq;

    // fill in the gapped positions. this is done for all nodes, including the root, the composition of non-root nodes is already correct
    for r in gaps.iter() {
      sequence[r.0..r.1].fill(alphabet.gap());
    }
    for r in unknown.iter() {
      // composition is already adjusted
      sequence[r.0..r.1].fill(alphabet.unknown());
    }
    if node.is_root {
      // if the node is the root, the composition is calculated from the full sequence
      *composition = Composition::with_sequence(sequence.iter().copied(), alphabet.chars(), alphabet.gap());
    }

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

    if !node.is_root {
      seq.sequence = seq![];
    }
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
  mut visitor: impl FnMut(&GraphNodeForward<NodeAncestral, EdgeAncestral, ()>, &Seq),
) -> Result<(), Report> {
  graph.iter_depth_first_preorder_forward(|node| {
    run_fitch_reconstruction(include_leaves, partitions, &mut visitor, &node);
  });
  Ok(())
}

fn run_fitch_reconstruction(
  include_leaves: bool,
  partitions: &[Arc<RwLock<PartitionFitch>>],
  mut visitor: impl FnMut(&GraphNodeForward<NodeAncestral, EdgeAncestral, ()>, &Seq),
  node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
) -> bool {
  if !include_leaves && node.is_leaf {
    return true;
  }

  for partition in partitions {
    let alphabet = &partition.read_arc().alphabet.clone(); // TODO: avoid clone

    let mut sequence = if !node.is_root {
      let partition = partition.read_arc();
      let (parent, edge) = get_exactly_one(&node.parent_keys).unwrap();
      let mut sequence = partition.nodes[parent].seq.sequence.clone();
      let edge_part = &partition.edges[edge];

      for sub in &edge_part.subs {
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

    visitor(node, &seq.sequence);
  }
  false
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
