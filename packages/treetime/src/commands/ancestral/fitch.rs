#![allow(dead_code)]

use crate::alphabet::alphabet::{FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward};
use crate::io::fasta::FastaRecord;
use crate::representation::graph_ancestral::{EdgeAncestral, GraphAncestral, NodeAncestral};
use crate::representation::graph_sparse::{
  Deletion, MarginalSparseSeqDistribution, ParsimonySeqDistribution, SparseEdgePartition, SparseNodePartition,
  SparseSeqInfo,
};
use crate::representation::partition_compressed::PartitionCompressed;
use crate::representation::partition_parsimony::PartitionParsimonyNew;
use crate::representation::seq::Seq;
use crate::representation::state_set::BitSet128;
use crate::representation::state_set::{StateSet, StateSetStatus};
use crate::seq::composition::Composition;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use crate::utils::container::get_exactly_one;
use crate::utils::interval::range::range_contains;
use crate::utils::interval::range_complement::range_complement;
use crate::utils::interval::range_difference::range_difference;
use crate::utils::interval::range_intersection::{range_intersection, range_intersection_iter};
use crate::{make_error, make_report, seq, stateset};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use maplit::btreemap;
use parking_lot::RwLock;
use std::sync::Arc;

fn attach_seqs_to_graph<P: PartitionCompressed>(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<P>>],
  aln: &[FastaRecord],
) -> Result<(), Report> {
  for leaf in graph.get_leaves() {
    let leaf_key = leaf.read_arc().key();
    let mut leaf = leaf.read_arc().payload().write_arc();

    let leaf_name = leaf.name.as_ref().ok_or_else(|| {
      make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
    })?.to_owned();

    let leaf_fasta = aln
      .iter()
      .find(|fasta| fasta.seq_name == leaf_name)
      // TODO: we could optionally emit a warning here and continue with a sequence that is missing
      .ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;

    *leaf = NodeAncestral {
      name: Some(leaf_name),
      desc: leaf_fasta.desc.clone(),
    };

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

fn fitch_backward<P: PartitionCompressed>(graph: &GraphAncestral, partitions: &[Arc<RwLock<P>>]) {
  graph.par_iter_breadth_first_backward(|node| {
    run_fitch_backward(partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
}

fn run_fitch_backward<P: PartitionCompressed>(
  partitions: &[Arc<RwLock<P>>],
  node: &GraphNodeBackward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  if node.is_leaf {
    return Ok(());
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
    let mut seq_dis = ParsimonySeqDistribution {
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

  Ok(())
}

fn fitch_forward<P: PartitionCompressed>(graph: &GraphAncestral, partitions: &[Arc<RwLock<P>>]) {
  graph.par_iter_breadth_first_forward(|node| {
    run_fitch_forward(partitions, &node).unwrap();
    GraphTraversalContinuation::Continue
  });
}

fn run_fitch_forward<P: PartitionCompressed>(
  partitions: &[Arc<RwLock<P>>],
  node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
) -> Result<(), Report> {
  for partition in partitions {
    let mut partition = partition.write_arc();
    let alphabet = &partition.alphabet().clone(); // TODO: avoid clone

    let mut node_data = partition.nodes_mut().remove(&node.key).unwrap();

    if node.is_root {
      let SparseSeqInfo {
        gaps,
        sequence,
        fitch: ParsimonySeqDistribution {
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
        fitch: ParsimonySeqDistribution {
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

fn fitch_cleanup<P: PartitionCompressed>(graph: &GraphAncestral, partitions: &[Arc<RwLock<P>>]) {
  graph.par_iter_breadth_first_forward(|node| run_fitch_forward_cleanup(&node, partitions));
}

fn run_fitch_forward_cleanup<P: PartitionCompressed>(
  node: &GraphNodeForward<NodeAncestral, EdgeAncestral, ()>,
  partitions: &[Arc<RwLock<P>>],
) -> GraphTraversalContinuation {
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

pub fn compress_sequences<P: PartitionCompressed>(
  graph: &GraphAncestral,
  partitions: &[Arc<RwLock<P>>],
  aln: &[FastaRecord],
) -> Result<(), Report> {
  attach_seqs_to_graph(graph, partitions, aln)?;
  fitch_backward(graph, partitions);
  fitch_forward(graph, partitions);
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
  partitions: &[Arc<RwLock<PartitionParsimonyNew>>],
  mut visitor: impl FnMut(&GraphNodeForward<NodeAncestral, EdgeAncestral, ()>, &Seq),
) -> Result<(), Report> {
  graph.iter_depth_first_preorder_forward(|node| {
    run_fitch_reconstruction(include_leaves, partitions, &mut visitor, &node);
  });
  Ok(())
}

fn run_fitch_reconstruction(
  include_leaves: bool,
  partitions: &[Arc<RwLock<PartitionParsimonyNew>>],
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

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::Alphabet;
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::json::{JsonPretty, json_write_str};
  use crate::io::nwk::nwk_read_str;
  use crate::representation::graph_ancestral::GraphAncestral;
  use crate::representation::partition_parsimony::PartitionParsimonyNew;
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use maplit::btreemap;
  use parking_lot::RwLock;
  use std::collections::BTreeMap;
  use std::sync::Arc;

  lazy_static! {
    static ref NUC_ALPHABET: Alphabet = Alphabet::default();
  }

  #[test]
  fn test_ancestral_reconstruction_fitch() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#},
      &NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      ACAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >CD
      CCGGCCATGTATTG--
    "#},
      &NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions_parsimony = [Arc::new(RwLock::new(PartitionParsimonyNew {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions_parsimony, &aln)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_fitch(&graph, false, &partitions_parsimony, |node, seq| {
      actual.insert(node.payload.name.clone(), seq.to_string());
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[test]
  fn test_ancestral_reconstruction_fitch_with_leaves() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
      indoc! {r#"
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#},
      &NUC_ALPHABET,
    )?;

    let expected = read_many_fasta_str(
      indoc! {r#"
      >root
      ACAGCCATGTATTG--
      >AB
      ACATCCCTGTA-TG--
      >A
      ACATCGCCNNA--GAC
      >B
      GCATCCCTGTA-NG--
      >CD
      CCGGCCATGTATTG--
      >C
      CCGGCGATGTRTTG--
      >D
      TCGGCCGTGTRTTG--
    "#},
      &NUC_ALPHABET,
    )?
    .into_iter()
    .map(|fasta| (fasta.seq_name, fasta.seq))
    .collect::<BTreeMap<_, _>>();

    let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();

    let partitions_parsimony = [Arc::new(RwLock::new(PartitionParsimonyNew {
      index: 0,
      alphabet,
      length: get_common_length(&aln)?,
      nodes: btreemap! {},
      edges: btreemap! {},
    }))];

    compress_sequences(&graph, &partitions_parsimony, &aln)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_fitch(&graph, true, &partitions_parsimony, |node, seq| {
      actual.insert(node.payload.name.clone(), seq.to_string());
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  // #[test]
  // fn test_fitch_internals() -> Result<(), Report> {
  //   rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

  //   let aln = read_many_fasta_str(
  //     indoc! {r#"
  //     >root
  //     ACAGCCATGTATTG--
  //     >AB
  //     ACATCCCTGTA-TG--
  //     >A
  //     ACATCGCCNNA--GAC
  //     >B
  //     GCATCCCTGTA-NG--
  //     >CD
  //     CCGGCCATGTATTG--
  //     >C
  //     CCGGCGATGTRTTG--
  //     >D
  //     TCGGCCGTGTRTTG--
  //   "#},
  //     &NUC_ALPHABET,
  //   )?;

  //   let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

  //   let alphabet = Alphabet::default();
  //   let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln.clone())?];

  //   attach_seqs_to_graph(&graph, &partitions)?;

  //   let partitions = partitions.into_iter().map(PartitionParsimony::from).collect_vec();

  //   fitch_backward(&graph, &partitions);

  //   {
  //     let seq_info = util::get_root_seq_info(&graph);
  //     assert_eq!(
  //       vec![0, 2, 3, 5, 6],
  //       seq_info.fitch.variable.keys().copied().collect_vec()
  //     );
  //     assert_eq!(&"~C~~C~~TGTATTGAC", &seq_info.sequence.as_str());
  //   }

  //   fitch_forward(&graph, &partitions);

  //   {
  //     let seq_info = util::get_root_seq_info(&graph);
  //     assert_eq!(&aln[0].seq, &seq_info.sequence);
  //   }

  //   {
  //     let actual_muts = util::collect_muts(&graph);
  //     let expected_muts = btreemap! {
  //       "AB->A"    => vec!["C6G", "T8C"],
  //       "AB->B"    => vec!["A1G"],
  //       "root->AB" => vec!["G4T", "A7C"],
  //       "CD->C"    => vec!["C6G"],
  //       "CD->D"    => vec!["C1T", "A7G"],
  //       "root->CD" => vec!["A1C", "A3G"],
  //     };
  //     assert_eq!(
  //       json_write_str(&expected_muts, JsonPretty(true))?,
  //       json_write_str(&actual_muts, JsonPretty(true))?
  //     );
  //   }

  //   {
  //     let actual_indels = util::collect_indels(&graph);
  //     let expected_indels = btreemap! {
  //       "AB->A"     => vec!["12--13: T -> -", "14--16: -- -> AC"],
  //       "AB->B"     => vec![],
  //       "CD->C"     => vec![],
  //       "root->AB"  => vec!["11--12: T -> -"],
  //       "CD->D"     => vec![],
  //       "root->CD"  => vec![],
  //     };
  //     assert_eq!(
  //       json_write_str(&expected_indels, JsonPretty(true))?,
  //       json_write_str(&actual_indels, JsonPretty(true))?
  //     );
  //   }

  //   Ok(())
  // }

  // #[test]
  // fn test_fitch_complex_gaps() -> Result<(), Report> {
  //   rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
  //
  //   // test the cases where: a) deletions overlap, b) the root has a deletion, c) an inserted sequence is variable
  //   // in DE, position 3 is inserted, but it varies in D and E
  //   let aln = read_many_fasta_str(
  //     indoc! {r#"
  //     >root
  //     TC-AG
  //     >AB
  //     TC-AG
  //     >A
  //     NC--G
  //     >B
  //     T--AG
  //     >CDE
  //     TG-TG
  //     >C
  //     TR-TG
  //     >DE
  //     TGTTG
  //     >D
  //     TGTTG
  //     >E
  //     TGCCG
  //   "#},
  //     &NUC_ALPHABET,
  //   )?;
  //
  //   let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,(D:0.05,E:0.03)DE:0.01)CDE:0.05)root:0.01;")?;
  //
  //   let alphabet = Alphabet::default();
  //   let partitions = vec![PartitionParsimonyWithAln::new(alphabet.clone(), aln.clone())?];
  //
  //   attach_seqs_to_graph(&graph, &partitions)?;
  //
  //   let partitions = partitions.into_iter().map(PartitionParsimony::from).collect_vec();
  //
  //   fitch_backward(&graph, &partitions);
  //
  //   fitch_forward(&graph, &partitions);
  //
  //   {
  //     let seq_info = util::get_root_seq_info(&graph);
  //     assert_eq!(&aln[0].seq, &seq_info.sequence);
  //   }
  //
  //   {
  //     let actual_muts = util::collect_muts(&graph);
  //     let expected_muts = btreemap! {
  //       "AB->A"     => vec![],
  //       "AB->B"     => vec![],
  //       "root->AB"  => vec![],
  //       "CDE->C"    => vec![],
  //       "CDE->DE"   => vec![],
  //       "DE->D"     => vec!["C3T"],
  //       "DE->E"     => vec!["T4C"],
  //       "root->CDE" => vec!["C2G", "A4T"],
  //     };
  //     assert_eq!(
  //       json_write_str(&expected_muts, JsonPretty(true))?,
  //       json_write_str(&actual_muts, JsonPretty(true))?
  //     );
  //   }
  //
  //   {
  //     let actual_indels = util::collect_indels(&graph);
  //     let expected_indels = btreemap! {
  //       "AB->A"     => vec!["3--4: A -> -"],
  //       "AB->B"     => vec!["1--2: C -> -"],
  //       "root->AB"  => vec![],
  //       "CDE->C"    => vec![],
  //       "CDE->DE"   => vec!["2--3: - -> C"],
  //       "DE->D"     => vec![],
  //       "DE->E"     => vec![],
  //       "root->CDE" => vec![],
  //     };
  //     assert_eq!(
  //       json_write_str(&expected_indels, JsonPretty(true))?,
  //       json_write_str(&actual_indels, JsonPretty(true))?
  //     );
  //   }
  //
  //   for node in graph.get_nodes() {
  //     let seq_info = util::get_root_seq_info(&graph);
  //     let composition = Composition::with_sequence(seq_info.sequence.iter().copied(), alphabet.chars(), alphabet.gap());
  //     assert_eq!(&seq_info.composition, &composition);
  //   }
  //
  //   Ok(())
  // }
  //
  // #[test]
  // fn test_fitch_polytomy() -> Result<(), Report> {
  //   rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
  //
  //   // test the cases where: a) deletions overlap, b) the root has a deletion, c) an inserted sequence is variable
  //   // in DE, position 3 is inserted, but it varies in D and E
  //   let aln = read_many_fasta_str(
  //     indoc! {r#"
  //     >root
  //     TC-AG
  //     >AB
  //     TC-AG
  //     >A
  //     NC--G
  //     >B
  //     T--AG
  //     >CDE
  //     TG-CG
  //     >C
  //     TR-TG
  //     >D
  //     TGTTG
  //     >E
  //     TGCCG
  //   "#},
  //     &NUC_ALPHABET,
  //   )?;
  //
  //   let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.05,E:0.03)CDE:0.05)root:0.01;")?;
  //
  //   let alphabet = Alphabet::default();
  //   let partitions = vec![PartitionParsimonyWithAln::new(alphabet.clone(), aln.clone())?];
  //
  //   attach_seqs_to_graph(&graph, &partitions)?;
  //
  //   let partitions = partitions.into_iter().map(PartitionParsimony::from).collect_vec();
  //
  //   fitch_backward(&graph, &partitions);
  //
  //   fitch_forward(&graph, &partitions);
  //
  //   {
  //     let seq_info = util::get_root_seq_info(&graph);
  //     assert_eq!(&aln[0].seq, &seq_info.sequence);
  //   }
  //
  //   {
  //     let actual_muts = util::collect_muts(&graph);
  //     let expected_muts = btreemap! {
  //       "AB->A"     => vec![],
  //       "AB->B"     => vec![],
  //       "root->AB"  => vec![],
  //       "CDE->C"    => vec!["C4T"],
  //       "CDE->D"    => vec!["C3T", "C4T"],
  //       "CDE->E"    => vec![],
  //       "root->CDE" => vec!["A4C","C2G"],
  //     };
  //     assert_eq!(
  //       json_write_str(&expected_muts, JsonPretty(true))?,
  //       json_write_str(&actual_muts, JsonPretty(true))?
  //     );
  //   }
  //
  //   {
  //     let actual_indels = util::collect_indels(&graph);
  //     let expected_indels = btreemap! {
  //       "AB->A"     => vec!["3--4: A -> -"],
  //       "AB->B"     => vec!["1--2: C -> -"],
  //       "root->AB"  => vec![],
  //       "CDE->C"    => vec!["2--3: C -> -"],
  //       "CDE->D"    => vec![],
  //       "CDE->E"    => vec![],
  //       "root->CDE" => vec!["2--3: - -> C"],
  //     };
  //     assert_eq!(
  //       json_write_str(&expected_indels, JsonPretty(true))?,
  //       json_write_str(&actual_indels, JsonPretty(true))?
  //     );
  //   }
  //
  //   for node in graph.get_nodes() {
  //     let seq_info = &node.read_arc().payload().read_arc().sparse_partitions[0].seq;
  //     let composition = Composition::with_sequence(seq_info.sequence.iter().copied(), alphabet.chars(), alphabet.gap());
  //     assert_eq!(&seq_info.composition, &composition);
  //   }
  //
  //   Ok(())
  // }
  //
  // #[test]
  // fn test_fitch_diverse_data() -> Result<(), Report> {
  //   rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
  //
  //   let alphabet = Alphabet::default();
  //   let aln = read_many_fasta(&["../../data/lassa/L/50/aln.fasta.xz"], &alphabet)?;
  //
  //   let graph: GraphAncestral = nwk_read_file("../../data/lassa/L/50/tree.nwk")?;
  //   let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln.clone())?];
  //   let partitions = compress_sequences(&graph, partitions)?;
  //
  //   // perform the ancestral reconstruction and save the sequences in a map
  //   let mut rec_seq = BTreeMap::new();
  //   ancestral_reconstruction_fitch(&graph, true, &partitions, |node, seq| {
  //     rec_seq.insert(node.name.clone().unwrap(), seq.to_string());
  //   })?;
  //
  //   let alphabet = Alphabet::default();
  //   for node in graph.get_inner_nodes() {
  //     let seq_info = &node.read_arc().payload().read_arc().sparse_partitions[0].seq;
  //     let node_name = node.read_arc().payload().read_arc().name.clone().unwrap();
  //     let composition = Composition::with_sequence(rec_seq[&node_name].chars(), alphabet.chars(), alphabet.gap());
  //     assert_eq!(&seq_info.composition, &composition);
  //   }
  //
  //   let alphabet = Alphabet::default();
  //   for node in graph.get_leaves() {
  //     let seq_info = &node.read_arc().payload().read_arc().sparse_partitions[0].seq;
  //     let node_name = node.read_arc().payload().read_arc().name.clone().unwrap();
  //     let input_seq = aln
  //       .iter()
  //       .find(|record| record.seq_name == node_name)
  //       .unwrap()
  //       .seq
  //       .iter()
  //       .copied();
  //     let input_seq_str: String = input_seq.clone().map(|c| c.to_string()).collect();
  //     assert_eq!(&input_seq_str, &rec_seq[&node_name]);
  //   }
  //   Ok(())
  // }
  //
  // mod util {
  //   use super::*;
  //
  //   pub fn get_root_seq_info(graph: &GraphAncestral) -> SparseSeqInfo {
  //     graph
  //       .get_exactly_one_root()
  //       .unwrap()
  //       .read_arc()
  //       .payload()
  //       .read_arc()
  //       .sparse_partitions[0]
  //       .seq
  //       .clone()
  //   }
  //
  //   pub fn collect_muts(graph: &GraphAncestral) -> BTreeMap<String, Vec<String>> {
  //     let actual_muts: BTreeMap<_, _> = graph
  //       .get_edges()
  //       .iter()
  //       .enumerate()
  //       .map(|(i, e)| {
  //         let get_name = |node_id| {
  //           let node = graph.get_node(node_id).unwrap().read_arc();
  //           node.payload().read_arc().name().unwrap().as_ref().to_owned()
  //         };
  //
  //         let src = get_name(e.read_arc().source());
  //         let tar = get_name(e.read_arc().target());
  //
  //         (
  //           format!("{src}->{tar}"),
  //           e.read_arc().payload().read_arc().sparse_partitions[0]
  //             .subs
  //             .iter()
  //             .map(ToString::to_string)
  //             .collect_vec(),
  //         )
  //       })
  //       .collect();
  //     actual_muts
  //   }
  //
  //   pub fn collect_indels(graph: &GraphAncestral) -> BTreeMap<String, Vec<String>> {
  //     let actual_indels: BTreeMap<_, _> = graph
  //       .get_edges()
  //       .iter()
  //       .map(|e| {
  //         let get_name = |node_id| {
  //           let node = graph.get_node(node_id).unwrap().read_arc();
  //           node.payload().read_arc().name().unwrap().as_ref().to_owned()
  //         };
  //
  //         let src = get_name(e.read_arc().source());
  //         let tar = get_name(e.read_arc().target());
  //
  //         (
  //           format!("{src}->{tar}"),
  //           e.read_arc().payload().read_arc().sparse_partitions[0]
  //             .indels
  //             .iter()
  //             .map(ToString::to_string)
  //             .collect_vec(),
  //         )
  //       })
  //       .collect();
  //     actual_indels
  //   }
  // }
}
