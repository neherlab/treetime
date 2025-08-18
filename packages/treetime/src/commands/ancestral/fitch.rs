#![allow(dead_code)]
use crate::alphabet::alphabet::{FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::node::Named;
use crate::io::fasta::FastaRecord;
use crate::representation::edge_partition::EdgePartition;
use crate::representation::graph_sparse::{
  Deletion, ParsimonySeqDis, SparseNodePartition, SparseSeqDis, SparseSeqInfo,
};
use crate::representation::node_partition::NodePartition;
use crate::representation::partitions_parsimony::{PartitionParsimony, PartitionParsimonyWithAln};
use crate::representation::repr_graph::{ReprGraph, ReprNode};
use crate::representation::seq::Seq;
use crate::representation::state_set::BitSet128;
use crate::representation::state_set::{StateSet, StateSetStatus};
use crate::seq::composition::Composition;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use crate::utils::interval::range::range_contains;
use crate::utils::interval::range_complement::range_complement;
use crate::utils::interval::range_difference::range_difference;
use crate::utils::interval::range_intersection::{range_intersection, range_intersection_iter};
use crate::{make_error, make_internal_report, make_report, seq, stateset};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use maplit::btreemap;
use ndarray::AssignElem;

fn attach_seqs_to_graph(graph: &ReprGraph, partitions: &[PartitionParsimonyWithAln]) -> Result<(), Report> {
  for leaf in graph.get_leaves() {
    let mut leaf = leaf.read_arc().payload().write_arc();

    let leaf_name = leaf.get_name().wrap_err(
      "Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name."
    )?.to_owned();

    // FIXME: all descs are the same for fasta partitions, so the mutable assignment here is needlessly complicated
    let mut desc = None;

    let partitions = partitions
      .iter()
      .map(
        |PartitionParsimonyWithAln { alphabet, aln, length }| -> Result<NodePartition, Report> {
          // TODO(perf): this might be slow if there are many sequences
          let leaf_fasta = aln
            .iter()
            .find(|fasta| fasta.seq_name == leaf_name)
            .ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;

          desc.assign_elem(leaf_fasta.desc.clone());

          //TODO: we could optionally emit a warning here and continue with a sequence that is entire missing...

          NodePartition::sparse(&leaf_fasta.seq, alphabet)
        },
      )
      .try_collect()?;

    *leaf = ReprNode::new(Some(leaf_name), desc, partitions);
  }

  Ok(())
}

fn fitch_backwards(graph: &ReprGraph, sparse_partitions: &[PartitionParsimony]) {
  let n_partitions = sparse_partitions.len();

  graph.par_iter_breadth_first_backward(|mut node| {
    if node.is_leaf {
      return GraphTraversalContinuation::Continue;
    }

    // initialize sparse partition data structure on incoming edges (could be moved to children, writing to outgoing edges)
    for (_, edge) in &node.children {
      edge.write_arc().set_partitions(
        (0..n_partitions)
          .map(|_| EdgePartition::sparse().unwrap())
          .collect_vec(),
      );
    }

    for si in 0..n_partitions {
      // the following are convenience datastructures
      let PartitionParsimony { alphabet, length } = &sparse_partitions[si];

      let children = node
        .children
        .iter()
        .map(|(c, e)| (c.read_arc(), e.read_arc()))
        .collect_vec();

      let children = children
        .iter()
        .map(|(c, e)| {
          (
            &c.partition_at(si).as_sparse().unwrap().seq,
            e.partition_at(si).as_sparse().unwrap(),
          )
        })
        .collect_vec();

      let n_children = children.len();

      // Initialization of target data structure (could be done later)
      let mut seq_dis = ParsimonySeqDis {
        variable: btreemap! {},
        variable_indel: btreemap! {},
        composition: Composition::new(alphabet.chars(), alphabet.gap()),
      };

      // determine parts of the sequence that are unknown, gaps in all children
      let mut gaps = range_intersection_iter(children.iter().map(|(c, _)| &c.gaps)).collect_vec();
      let unknown = range_intersection_iter(children.iter().map(|(c, _)| &c.unknown)).collect_vec();
      // non_char are ranges that are either unknown or gaps in all children (note that can be different from the union of gaps and unknown)
      let non_char = range_intersection_iter(children.iter().map(|(c, _)| &c.non_char)).collect_vec();
      // calculate the complement of gaps for later look-up
      let non_gap = range_complement(&[(0, *length)], &[gaps.clone()]); // FIXME(perf): unnecessary clone

      // what follows could be a function that returns `sequence` and `variable`, takes as arguments children, non_char, alphabet
      let mut sequence = seq![FILL_CHAR; *length];
      for r in &non_char {
        sequence.apply_unknowns(*r, NON_CHAR);
      }

      // Process all positions where the children are variable.
      // Need to account for parts of the sequence transmitted along edges.
      let variable_positions = children
        .iter()
        .flat_map(|(c, _)| c.fitch.variable.keys().copied())
        .unique()
        .collect_vec();

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
            sequence.set_char(pos, state);
          },
          StateSetStatus::Ambiguous(_) => {
            // more than one possible states
            seq_dis.variable.insert(pos, intersection);
            sequence.set_char(pos, VARIABLE_CHAR);
          },
          StateSetStatus::Empty => {
            let union = StateSet::from_union(&child_profiles);
            seq_dis.variable.insert(pos, union);
            sequence.set_char(pos, VARIABLE_CHAR);
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
          if alphabet.is_canonical(child_state) {
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

      node.payload.set_partition_at(
        si,
        NodePartition::Sparse(SparseNodePartition {
          seq: SparseSeqInfo {
            gaps,
            unknown,
            non_char,
            fitch: seq_dis,
            sequence,
            composition: Composition::new(alphabet.chars(), alphabet.gap()),
          },
          profile: SparseSeqDis {
            variable: btreemap! {},
            variable_indel: btreemap! {},
            fixed: btreemap! {},
            fixed_counts: Composition::new(alphabet.chars(), alphabet.gap()),
            log_lh: 0.0,
          },
        }),
      );
    }

    GraphTraversalContinuation::Continue
  });
}

fn fitch_forward(graph: &ReprGraph, sparse_partitions: &[PartitionParsimony]) {
  let n_partitions = sparse_partitions.len();

  graph.par_iter_breadth_first_forward(|mut node| {
    for si in 0..n_partitions {
      let PartitionParsimony { alphabet, .. } = &sparse_partitions[si];
      let SparseSeqInfo {
        gaps,
        unknown,
        sequence,
        composition,
        non_char,
        fitch: ParsimonySeqDis {
          variable,
          variable_indel,
          ..
        },
      } = &mut node.payload.partition_at_mut(si).as_sparse_mut().unwrap().seq;

      if node.is_root {
        for (pos, states) in variable {
          sequence.set_char(*pos, states.get_one());
        }
        // process indels as majority rule at the root
        for (r, indel) in variable_indel.iter() {
          if indel.deleted > indel.present {
            gaps.push(*r);
          }
        }
      } else {
        let (parent, edge) = node
          .parents
          .first()
          .ok_or_else(|| make_internal_report!("Graphs with multiple parents per node are not yet supported"))
          .unwrap();

        let parent = parent.read_arc();
        let parent = &parent.partition_at(si).as_sparse().unwrap().seq;
        let mut edge = edge.write_arc();
        let edge = &mut edge.partition_at_mut(si).as_sparse_mut().unwrap();

        *composition = parent.composition.clone();

        // fill in the indeterminate positions by copying the parent (note that new gaps in the node will be introduced later)
        for r in non_char {
          sequence[r.0..r.1].clone_from_slice(&parent.sequence[r.0..r.1]);
        }

        // the following two loops modify the sequence, composition, and edge in place and process variable position.
        // for each variable position, pick a state or a mutation
        for (pos, states) in variable.iter_mut() {
          let parent_char = parent.sequence[*pos];
          if alphabet.is_canonical(parent_char) {
            if states.contains(parent_char) {
              sequence.set_char(*pos, parent_char);
            } else {
              let child_char = states.get_one();
              sequence.set_char(*pos, child_char);
              let m = Sub::new(parent_char, *pos, child_char).unwrap();
              m.check_determined(alphabet).unwrap();
              composition.add_sub(&m);
              edge.subs.push(m);
            }
          } else if alphabet.is_gap(parent_char) && !range_contains(gaps, *pos) {
            // if parent is gap, but child isn't, we need to resolve variable states
            sequence.set_char(*pos, states.get_one());
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
            let m = Sub::new(parent.sequence[pos], pos, node_nuc).unwrap();
            m.check_determined(alphabet).unwrap();
            composition.add_sub(&m);
            edge.subs.push(m);
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
              edge.indels.push(indel);
            }
          } else if gap_in_parent > 0 {
            // Add insertion if gap is present in parent.
            let indel = InDel::ins(*r, &sequence[r.0..r.1]);
            composition.add_indel(&indel);
            edge.indels.push(indel);
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
          edge.indels.push(indel);
        }

        // Process gaps in the parent that are not in the node (insertions)
        for r in range_difference(&parent.gaps, gaps) {
          if variable_indel.contains_key(&r) {
            // all gaps in variable_indel are already processed
            continue;
          }
          let indel = InDel::ins(r, &sequence[r.0..r.1]);
          composition.add_indel(&indel);
          edge.indels.push(indel);
        }
        for r in unknown.iter() {
          // this might result in compensating addition/deletions of Ns already present in the parent
          for pos in r.0..r.1 {
            composition.adjust_count(sequence[pos], -1);
          }
          composition.adjust_count(alphabet.unknown(), r.1 as isize - r.0 as isize);
        }
      }
      // fill in the gapped positions. this is done for all nodes, including the root, the composition of non-root nodes is already correct
      for r in gaps.iter() {
        sequence.apply_del(*r, alphabet.gap());
      }
      for r in unknown.iter() {
        // composition is already adjusted
        sequence.apply_unknowns(*r, alphabet.unknown());
      }
      if node.is_root {
        // if the node is the root, the composition is calculated from the full sequence
        *composition = Composition::with_sequence(sequence.iter().copied(), alphabet.chars(), alphabet.gap());
      }
    }
    GraphTraversalContinuation::Continue
  });
}

fn fitch_cleanup(graph: &ReprGraph) {
  graph.par_iter_breadth_first_forward(|mut node| {
    for part in node.payload.partitions_mut() {
      let seq = &mut part.as_sparse_mut().unwrap().seq;

      // delete the variable position everywhere instead of leaves
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
  });
}

pub fn compress_sequences(
  graph: &ReprGraph,
  partitions: Vec<PartitionParsimonyWithAln>,
) -> Result<Vec<PartitionParsimony>, Report> {
  attach_seqs_to_graph(graph, &partitions)?;

  let partitions = partitions.into_iter().map(PartitionParsimony::from).collect_vec();

  fitch_backwards(graph, &partitions);
  fitch_forward(graph, &partitions);
  fitch_cleanup(graph);

  Ok(partitions)
}

/// Reconstruct ancestral sequences using Fitch parsimony.
///
/// Calls visitor function for every ancestral node, providing the node itself and its reconstructed sequence.
/// Optionally reconstructs leaf sequences.
pub fn ancestral_reconstruction_fitch(
  graph: &ReprGraph,
  include_leaves: bool,
  partitions: &[PartitionParsimony],
  mut visitor: impl FnMut(&ReprNode, &Seq),
) -> Result<(), Report> {
  let n_partitions = partitions.len();

  graph.iter_depth_first_preorder_forward(|mut node| {
    if !include_leaves && node.is_leaf {
      return;
    }

    for si in 0..n_partitions {
      let PartitionParsimony { alphabet, .. } = &partitions[si];

      if !node.is_root {
        let (parent, edge) = node.get_exactly_one_parent().unwrap();

        let parent = parent.read_arc();
        let parent_sequence = &parent.partition_at(si).as_sparse().unwrap().seq.sequence;

        let edge = edge.read_arc();
        let edge_partition = &edge.partition_at(si).as_sparse().unwrap();
        let subs = &edge_partition.subs;
        let indels = &edge_partition.indels;

        let seq = &mut node.payload.partition_at_mut(si).as_sparse_mut().unwrap().seq;
        let sequence = &mut seq.sequence;
        *sequence = parent_sequence.clone();

        for sub in subs {
          sequence.apply_sub(sub);
        }
        for indel in indels {
          sequence.apply_indel(indel, alphabet.gap());
        }
      }

      {
        let seq = &mut node.payload.partition_at_mut(si).as_sparse_mut().unwrap().seq;
        let sequence = &mut seq.sequence;

        for r in &seq.unknown {
          sequence.apply_unknowns(*r, alphabet.unknown());
        }

        for (pos, states) in &seq.fitch.variable {
          sequence.set_char(*pos, alphabet.set_to_char(*states));
        }
      }

      let sequence = &node.payload.partition_at(si).as_sparse().unwrap().seq.sequence;
      visitor(&node.payload, &sequence);
    }
  });

  Ok(())
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
  use crate::graph::node::Named;
  use crate::io::fasta::{read_many_fasta, read_many_fasta_str};
  use crate::io::json::{JsonPretty, json_write_str};
  use crate::io::nwk::{nwk_read_file, nwk_read_str};
  use eyre::Report;
  use indoc::indoc;
  use lazy_static::lazy_static;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;

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

    let graph: ReprGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
    let partitions = compress_sequences(&graph, partitions)?;

    let mut actual = btreemap! {};
    ancestral_reconstruction_fitch(&graph, false, &partitions, |node, seq| {
      actual.insert(node.get_name().unwrap().to_owned(), seq.to_string());
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

    let graph: ReprGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
    let partitions = compress_sequences(&graph, partitions)?;

    let mut actual = btreemap! {};
    ancestral_reconstruction_fitch(&graph, true, &partitions, |node, seq| {
      actual.insert(node.get_name().unwrap().to_owned(), seq.to_string());
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  #[test]
  fn test_fitch_internals() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let aln = read_many_fasta_str(
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
    )?;

    let graph: ReprGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln.clone())?];

    attach_seqs_to_graph(&graph, &partitions)?;

    let partitions = partitions.into_iter().map(PartitionParsimony::from).collect_vec();

    fitch_backwards(&graph, &partitions);

    {
      let seq_info = util::get_root_seq_info(&graph);
      assert_eq!(
        vec![0, 2, 3, 5, 6],
        seq_info.fitch.variable.keys().copied().collect_vec()
      );
      assert_eq!(&"~C~~C~~TGTATTGAC", &seq_info.sequence.as_str());
    }

    fitch_forward(&graph, &partitions);

    {
      let seq_info = util::get_root_seq_info(&graph);
      assert_eq!(&aln[0].seq, &seq_info.sequence);
    }

    {
      let actual_muts = util::collect_muts(&graph);
      let expected_muts = btreemap! {
        "AB->A"    => vec!["C6G", "T8C"],
        "AB->B"    => vec!["A1G"],
        "root->AB" => vec!["G4T", "A7C"],
        "CD->C"    => vec!["C6G"],
        "CD->D"    => vec!["C1T", "A7G"],
        "root->CD" => vec!["A1C", "A3G"],
      };
      assert_eq!(
        json_write_str(&expected_muts, JsonPretty(true))?,
        json_write_str(&actual_muts, JsonPretty(true))?
      );
    }

    {
      let actual_indels = util::collect_indels(&graph);
      let expected_indels = btreemap! {
        "AB->A"     => vec!["12--13: T -> -", "14--16: -- -> AC"],
        "AB->B"     => vec![],
        "CD->C"     => vec![],
        "root->AB"  => vec!["11--12: T -> -"],
        "CD->D"     => vec![],
        "root->CD"  => vec![],
      };
      assert_eq!(
        json_write_str(&expected_indels, JsonPretty(true))?,
        json_write_str(&actual_indels, JsonPretty(true))?
      );
    }

    Ok(())
  }

  #[test]
  fn test_fitch_complex_gaps() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    // test the cases where: a) deletions overlap, b) the root has a deletion, c) an inserted sequence is variable
    // in DE, position 3 is inserted, but it varies in D and E
    let aln = read_many_fasta_str(
      indoc! {r#"
      >root
      TC-AG
      >AB
      TC-AG
      >A
      NC--G
      >B
      T--AG
      >CDE
      TG-TG
      >C
      TR-TG
      >DE
      TGTTG
      >D
      TGTTG
      >E
      TGCCG
    "#},
      &NUC_ALPHABET,
    )?;

    let graph: ReprGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,(D:0.05,E:0.03)DE:0.01)CDE:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet.clone(), aln.clone())?];

    attach_seqs_to_graph(&graph, &partitions)?;

    let partitions = partitions.into_iter().map(PartitionParsimony::from).collect_vec();

    fitch_backwards(&graph, &partitions);

    fitch_forward(&graph, &partitions);

    {
      let seq_info = util::get_root_seq_info(&graph);
      assert_eq!(&aln[0].seq, &seq_info.sequence);
    }

    {
      let actual_muts = util::collect_muts(&graph);
      let expected_muts = btreemap! {
        "AB->A"     => vec![],
        "AB->B"     => vec![],
        "root->AB"  => vec![],
        "CDE->C"    => vec![],
        "CDE->DE"   => vec![],
        "DE->D"     => vec!["C3T"],
        "DE->E"     => vec!["T4C"],
        "root->CDE" => vec!["C2G", "A4T"],
      };
      assert_eq!(
        json_write_str(&expected_muts, JsonPretty(true))?,
        json_write_str(&actual_muts, JsonPretty(true))?
      );
    }

    {
      let actual_indels = util::collect_indels(&graph);
      let expected_indels = btreemap! {
        "AB->A"     => vec!["3--4: A -> -"],
        "AB->B"     => vec!["1--2: C -> -"],
        "root->AB"  => vec![],
        "CDE->C"    => vec![],
        "CDE->DE"   => vec!["2--3: - -> C"],
        "DE->D"     => vec![],
        "DE->E"     => vec![],
        "root->CDE" => vec![],
      };
      assert_eq!(
        json_write_str(&expected_indels, JsonPretty(true))?,
        json_write_str(&actual_indels, JsonPretty(true))?
      );
    }

    for node in graph.get_nodes() {
      let seq_info = util::get_root_seq_info(&graph);
      let composition = Composition::with_sequence(seq_info.sequence.iter().copied(), alphabet.chars(), alphabet.gap());
      assert_eq!(&seq_info.composition, &composition);
    }

    Ok(())
  }

  #[test]
  fn test_fitch_polytomy() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    // test the cases where: a) deletions overlap, b) the root has a deletion, c) an inserted sequence is variable
    // in DE, position 3 is inserted, but it varies in D and E
    let aln = read_many_fasta_str(
      indoc! {r#"
      >root
      TC-AG
      >AB
      TC-AG
      >A
      NC--G
      >B
      T--AG
      >CDE
      TG-CG
      >C
      TR-TG
      >D
      TGTTG
      >E
      TGCCG
    "#},
      &NUC_ALPHABET,
    )?;

    let graph: ReprGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.05,E:0.03)CDE:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet.clone(), aln.clone())?];

    attach_seqs_to_graph(&graph, &partitions)?;

    let partitions = partitions.into_iter().map(PartitionParsimony::from).collect_vec();

    fitch_backwards(&graph, &partitions);

    fitch_forward(&graph, &partitions);

    {
      let seq_info = util::get_root_seq_info(&graph);
      assert_eq!(&aln[0].seq, &seq_info.sequence);
    }

    {
      let actual_muts = util::collect_muts(&graph);
      let expected_muts = btreemap! {
        "AB->A"     => vec![],
        "AB->B"     => vec![],
        "root->AB"  => vec![],
        "CDE->C"    => vec!["C4T"],
        "CDE->D"    => vec!["C3T", "C4T"],
        "CDE->E"    => vec![],
        "root->CDE" => vec!["A4C","C2G"],
      };
      assert_eq!(
        json_write_str(&expected_muts, JsonPretty(true))?,
        json_write_str(&actual_muts, JsonPretty(true))?
      );
    }

    {
      let actual_indels = util::collect_indels(&graph);
      let expected_indels = btreemap! {
        "AB->A"     => vec!["3--4: A -> -"],
        "AB->B"     => vec!["1--2: C -> -"],
        "root->AB"  => vec![],
        "CDE->C"    => vec!["2--3: C -> -"],
        "CDE->D"    => vec![],
        "CDE->E"    => vec![],
        "root->CDE" => vec!["2--3: - -> C"],
      };
      assert_eq!(
        json_write_str(&expected_indels, JsonPretty(true))?,
        json_write_str(&actual_indels, JsonPretty(true))?
      );
    }

    for node in graph.get_nodes() {
      let node = node.read_arc();
      let node = node.payload().read_arc();
      let seq_info = &node.partition_at(0).as_sparse()?.seq;
      let composition = Composition::with_sequence(seq_info.sequence.iter().copied(), alphabet.chars(), alphabet.gap());
      assert_eq!(&seq_info.composition, &composition);
    }

    Ok(())
  }

  #[test]
  fn test_fitch_diverse_data() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let alphabet = Alphabet::default();
    let aln = read_many_fasta(&["../../data/lassa/L/50/aln.fasta.xz"], &alphabet)?;

    let graph: ReprGraph = nwk_read_file("../../data/lassa/L/50/tree.nwk")?;
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln.clone())?];
    let partitions = compress_sequences(&graph, partitions)?;

    let mut reconstructed = btreemap! {};
    ancestral_reconstruction_fitch(&graph, true, &partitions, |node, seq| {
      reconstructed.insert(node.get_name().unwrap().to_owned(), seq.to_string());
    })?;

    let alphabet = Alphabet::default();
    for node in graph.get_inner_nodes() {
      let node = node.read_arc();
      let node = node.payload().read_arc();
      let seq_info = &node.partition_at(0).as_sparse()?.seq;
      let node_name = node.get_name()?;
      let composition = Composition::with_sequence(reconstructed[node_name].chars(), alphabet.chars(), alphabet.gap());
      assert_eq!(&seq_info.composition, &composition);
    }

    let alphabet = Alphabet::default();
    for node in graph.get_leaves() {
      let node = node.read_arc();
      let node = node.payload().read_arc();
      let seq_info = &node.partition_at(0).as_sparse()?.seq;
      let node_name = node.get_name()?;
      let input_seq = aln
        .iter()
        .find(|record| record.seq_name == node_name)
        .unwrap()
        .seq
        .iter()
        .copied();
      let input_seq_str: String = input_seq.clone().map(|c| c.to_string()).collect();
      assert_eq!(&input_seq_str, &reconstructed[node_name]);
    }
    Ok(())
  }

  mod util {
    use super::*;

    pub fn get_root_seq_info(graph: &ReprGraph) -> SparseSeqInfo {
      graph
        .get_exactly_one_root()
        .unwrap()
        .read_arc()
        .payload()
        .read_arc()
        .partition_at(0)
        .as_sparse()
        .unwrap()
        .seq
        .clone()
    }

    pub fn collect_muts(graph: &ReprGraph) -> BTreeMap<String, Vec<String>> {
      let actual_muts: BTreeMap<_, _> = graph
        .get_edges()
        .iter()
        .enumerate()
        .map(|(i, e)| {
          let get_name = |node_id| {
            let node = graph.get_node(node_id).unwrap().read_arc();
            node.payload().read_arc().get_name().unwrap().to_owned()
          };

          let src = get_name(e.read_arc().source());
          let tar = get_name(e.read_arc().target());

          (
            format!("{src}->{tar}"),
            e.read_arc()
              .payload()
              .read_arc()
              .partition_at(0)
              .as_sparse()
              .unwrap()
              .subs
              .iter()
              .map(ToString::to_string)
              .collect_vec(),
          )
        })
        .collect();
      actual_muts
    }

    pub fn collect_indels(graph: &ReprGraph) -> BTreeMap<String, Vec<String>> {
      let actual_indels: BTreeMap<_, _> = graph
        .get_edges()
        .iter()
        .map(|e| {
          let get_name = |node_id| {
            let node = graph.get_node(node_id).unwrap().read_arc();
            node.payload().read_arc().get_name().unwrap().to_owned()
          };

          let src = get_name(e.read_arc().source());
          let tar = get_name(e.read_arc().target());

          (
            format!("{src}->{tar}"),
            e.read_arc()
              .payload()
              .read_arc()
              .partition_at(0)
              .as_sparse()
              .unwrap()
              .indels
              .iter()
              .map(ToString::to_string)
              .collect_vec(),
          )
        })
        .collect();
      actual_indels
    }
  }
}
