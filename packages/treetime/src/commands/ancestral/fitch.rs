#![allow(dead_code)]
use crate::alphabet::alphabet::{FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::io::fasta::FastaRecord;
use crate::representation::graph_sparse::{
  Deletion, ParsimonySeqDis, SparseGraph, SparseNode, SparseSeqDis, SparseSeqEdge, SparseSeqInfo, SparseSeqNode,
};
use crate::representation::partitions_parsimony::{PartitionParsimony, PartitionParsimonyWithAln};
use crate::representation::state_set::BitSet128;
use crate::representation::state_set::{StateSet, StateSetStatus};
use crate::seq::composition::Composition;
use crate::seq::indel::InDel;
use crate::seq::mutation::Sub;
use crate::utils::interval::range::range_contains;
use crate::utils::interval::range_complement::range_complement;
use crate::utils::interval::range_difference::range_difference;
use crate::utils::interval::range_intersection::{range_intersection, range_intersection_iter};
use crate::{make_error, make_internal_report, make_report, stateset};
use eyre::{Report, WrapErr};
use itertools::Itertools;
use maplit::btreemap;
use ndarray::AssignElem;

fn attach_seqs_to_graph(graph: &SparseGraph, partitions: &[PartitionParsimonyWithAln]) -> Result<(), Report> {
  for leaf in graph.get_leaves() {
    let mut leaf = leaf.read_arc().payload().write_arc();

    let leaf_name = leaf.name.as_ref().ok_or_else(|| {
      make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
    })?.to_owned();

    // FIXME: all descs are the same for fasta partitions, so the mutable assignment here is needlessly complicated
    let mut desc = None;

    let sparse_partitions = partitions
      .iter()
      .map(|PartitionParsimonyWithAln { alphabet, aln, length }| {
        // TODO(perf): this might be slow if there are many sequences
        let leaf_fasta = aln
          .iter()
          .find(|fasta| fasta.seq_name == leaf_name)
          .ok_or_else(|| make_report!("Leaf sequence not found: '{leaf_name}'"))?;

        desc.assign_elem(leaf_fasta.desc.clone());

        //TODO: we could optionally emit a warning here and continue with a sequence that is entire missing...

        // TODO(perf): unnecessary copy of sequence data. Neither String, nor &[char] works well for us, it seems.
        // We probably want a custom class for sequences. Sequences should be instantiated in the fasta parser and
        // never need a copy like here.
        let sequence = leaf_fasta.seq.chars().collect::<Vec<_>>();

        SparseSeqNode::new(&sequence, alphabet)
      })
      .collect::<Result<_, Report>>()?;

    *leaf = SparseNode {
      name: Some(leaf_name),
      desc,
      sparse_partitions,
    }
  }

  Ok(())
}

fn fitch_backwards(graph: &SparseGraph, sparse_partitions: &[PartitionParsimony]) {
  let n_partitions = sparse_partitions.len();

  graph.par_iter_breadth_first_backward(|mut node| {
    if node.is_leaf {
      return GraphTraversalContinuation::Continue;
    }

    // initialize sparse partition data structure on incoming edges (could be moved to children, writing to outgoing edges)
    for (_, edge) in &node.children {
      edge.write_arc().sparse_partitions = vec![SparseSeqEdge::default(); n_partitions];
    }

    for si in 0..n_partitions {
      let PartitionParsimony { alphabet, length } = &sparse_partitions[si];

      let children = node
        .children
        .iter()
        .map(|(c, e)| (c.read_arc(), e.read_arc()))
        .collect_vec();

      let children = children
        .iter()
        .map(|(c, e)| (&c.sparse_partitions[si].seq, &e.sparse_partitions[si]))
        .collect_vec();

      let n_children = children.len();

      // determine parts of the sequence that are unknown, gaps in all children
      let mut gaps = range_intersection_iter(children.iter().map(|(c, _)| &c.gaps)).collect_vec();
      let unknown = range_intersection_iter(children.iter().map(|(c, _)| &c.unknown)).collect_vec();
      // non_char are ranges that are either unknown or gaps in all children (note that can be different from the union of gaps and unknown)
      let non_char = range_intersection_iter(children.iter().map(|(c, _)| &c.non_char)).collect_vec();
      // calculate the complement of gaps for later look-up
      let non_gap = range_complement(&[(0, *length)], &[gaps.clone()]); // FIXME(perf): unnecessary clone

      //
      let mut seq_dis = ParsimonySeqDis {
        variable: btreemap! {},
        variable_indel: btreemap! {},
        composition: Composition::new(alphabet.chars(), alphabet.gap()),
      };
      let mut sequence = vec![FILL_CHAR; *length];

      for r in &non_char {
        sequence[r.0..r.1].fill(NON_CHAR);
      }

      // Process all positions where the children are variable.
      // Need to account for parts of the sequence transmitted along edges.
      let variable_positions = children
        .iter()
        .flat_map(|(c, e)| c.fitch.variable.keys().copied())
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
            sequence[pos] = state;
          }
          StateSetStatus::Ambiguous(_) => {
            // more than one possible states
            seq_dis.variable.insert(pos, intersection);
            sequence[pos] = VARIABLE_CHAR;
          }
          StateSetStatus::Empty => {
            let union = StateSet::from_union(&child_profiles);
            seq_dis.variable.insert(pos, union);
            sequence[pos] = VARIABLE_CHAR;
          }
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
              *seq_dis.variable.entry(pos).or_insert(stateset! {*parent_state}) += child_state;
              *parent_state = VARIABLE_CHAR;
            }
          }
        }
      }

      // Process insertions and deletions.

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

      node.payload.sparse_partitions.push(SparseSeqNode {
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
      });
    }

    GraphTraversalContinuation::Continue
  });
}

fn fitch_forward(graph: &SparseGraph, sparse_partitions: &[PartitionParsimony]) {
  let n_partitions = sparse_partitions.len();

  graph.par_iter_breadth_first_forward(|mut node| {
    for si in 0..n_partitions {
      let PartitionParsimony { alphabet, .. } = &sparse_partitions[si];
      let SparseSeqInfo {
        gaps,
        sequence,
        composition,
        non_char,
        fitch: ParsimonySeqDis {
          variable,
          variable_indel,
          ..
        },
        ..
      } = &mut node.payload.sparse_partitions[si].seq;

      if node.is_root {
        for (pos, states) in variable {
          sequence[*pos] = states.get_one();
        }
        // process indels as majority rule at the root
        for (r, indel) in variable_indel.iter() {
          if indel.deleted > indel.present {
            gaps.push(*r);
          }
        }
        for r in gaps.iter() {
          sequence[r.0..r.1].fill(alphabet.gap());
        }
        *composition = Composition::with_sequence(sequence.iter().copied(), alphabet.chars(), alphabet.gap());
      } else {
        let (parent, edge) = node
          .parents
          .first()
          .ok_or_else(|| make_internal_report!("Graphs with multiple parents per node are not yet supported"))
          .unwrap();

        let parent = &parent.read_arc().sparse_partitions[si].seq;
        let edge = &mut edge.write_arc().sparse_partitions[si];

        *composition = parent.composition.clone();

        // fill in the indeterminate positions
        for r in non_char {
          sequence[r.0..r.1].clone_from_slice(&parent.sequence[r.0..r.1]);
        }

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
              let m = Sub {
                pos: *pos,
                qry: cnuc,
                reff: pnuc,
              };
              composition.add_sub(&m);
              edge.subs.push(m);
            }
          } else if alphabet.is_gap(pnuc) && !range_contains(gaps, *pos) {
            // if parent is gap, but child isn't, we need to resolve variable states
            sequence[*pos] = states.get_one();
          }
        }

        for &pos in parent.fitch.variable.keys() {
          if variable.contains_key(&pos) {
            continue;
          }

          // NOTE: access to full_seq would not be necessary if we had saved the
          // child state of variable positions in the backward pass
          let node_nuc = sequence[pos];
          if alphabet.is_canonical(node_nuc) && parent.sequence[pos] != node_nuc {
            let m = Sub {
              pos,
              qry: node_nuc,
              reff: parent.sequence[pos],
            };
            composition.add_sub(&m);
            edge.subs.push(m);
          }
        }

        // Process indels. Gaps where the children disagree, need to be decided by also looking at parent.
        for (r, indel) in variable_indel.iter() {
          let gap_in_parent = if parent.gaps.contains(r) { 1 } else { 0 };
          if indel.deleted + gap_in_parent > indel.present {
            gaps.push(*r);
            if gap_in_parent == 0 {
              // If the gap is not in parent, add deletion.
              // the sequence that is deleted is the sequence of the parent
              let indel = InDel {
                range: *r,
                seq: parent.sequence[r.0..r.1].to_owned(),
                deletion: true,
              };
              composition.add_indel(&indel);
              edge.indels.push(indel);
            }
          } else if gap_in_parent > 0 {
            // Add insertion if gap is present in parent.
            let indel = InDel {
              range: *r,
              seq: sequence[r.0..r.1].to_owned(),
              deletion: false,
            };
            composition.add_indel(&indel);
            edge.indels.push(indel);
          }
        }

        // Process consensus gaps in the node that are not in the parent (deletions)
        for r in range_difference(gaps, &parent.gaps) {
          let indel = InDel {
            range: r,
            seq: sequence[r.0..r.1].to_owned(),
            deletion: true,
          };
          composition.add_indel(&indel);
          edge.indels.push(indel);
        }

        // Process gaps in the parent that are not in the node (insertions)
        for r in range_difference(&parent.gaps, gaps) {
          let indel = InDel {
            range: r,
            seq: sequence[r.0..r.1].to_owned(),
            deletion: false,
          };
          composition.add_indel(&indel);
          edge.indels.push(indel);
        }
      }
    }
    GraphTraversalContinuation::Continue
  });
}

fn fitch_cleanup(graph: &SparseGraph) {
  graph.par_iter_breadth_first_forward(|mut node| {
    for SparseSeqNode { seq, .. } in &mut node.payload.sparse_partitions {
      // delete the variable position everywhere instead of leaves
      if !node.is_leaf {
        seq.fitch.variable = btreemap! {};
      }

      // remove the undetermined counts from the counts of fixed positions
      for r in &seq.unknown {
        for pos in r.0..r.1 {
          seq.composition.adjust_count(seq.sequence[pos], -1);
        }
      }

      seq.fitch.composition = seq.composition.clone();
      for p in seq.fitch.variable.values() {
        if let Some(state) = p.get_one_maybe() {
          seq.fitch.composition.adjust_count(state, -1);
        }
      }

      if !node.is_root {
        seq.sequence = vec![];
      }
    }

    GraphTraversalContinuation::Continue
  });
}

pub fn compress_sequences(
  graph: &SparseGraph,
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
  graph: &SparseGraph,
  include_leaves: bool,
  partitions: &[PartitionParsimony],
  mut visitor: impl FnMut(&SparseNode, &[char]),
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
        let parent_seq = &parent.read_arc().sparse_partitions[si].seq.sequence;
        let edge_part = &edge.read_arc().sparse_partitions[si];

        node.payload.sparse_partitions[si].seq.sequence = parent_seq.clone();

        for sub in &edge_part.subs {
          node.payload.sparse_partitions[si].seq.sequence[sub.pos] = sub.qry;
        }

        for indel in &edge_part.indels {
          if indel.deletion {
            node.payload.sparse_partitions[si].seq.sequence[indel.range.0..indel.range.1].fill(alphabet.gap());
          } else {
            node.payload.sparse_partitions[si].seq.sequence[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
          }
        }
      }

      let seq = &mut node.payload.sparse_partitions[si].seq;

      for r in &mut seq.unknown {
        seq.sequence[r.0..r.1].fill(alphabet.unknown());
      }

      for (pos, states) in &mut seq.fitch.variable {
        seq.sequence[*pos] = alphabet.set_to_char(*states);
      }

      visitor(&node.payload, &node.payload.sparse_partitions[si].seq.sequence);
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
    }
  }
  .wrap_err("When calculating length of sequences")
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::graph::node::Named;
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::json::{json_write_str, JsonPretty};
  use crate::io::nwk::nwk_read_str;
  use crate::utils::string::vec_to_string;
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

    let graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
    let partitions = compress_sequences(&graph, partitions)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_fitch(&graph, false, &partitions, |node, seq| {
      actual.insert(node.name.clone(), vec_to_string(seq.to_owned()));
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

    let graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln)?];
    let partitions = compress_sequences(&graph, partitions)?;

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_fitch(&graph, true, &partitions, |node, seq| {
      actual.insert(node.name.clone(), vec_to_string(seq.to_owned()));
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

    let graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let alphabet = Alphabet::default();
    let partitions = vec![PartitionParsimonyWithAln::new(alphabet, aln.clone())?];

    attach_seqs_to_graph(&graph, &partitions)?;

    let partitions = partitions.into_iter().map(PartitionParsimony::from).collect_vec();

    fitch_backwards(&graph, &partitions);

    {
      let seq_info = &graph
        .get_exactly_one_root()?
        .read_arc()
        .payload()
        .read_arc()
        .sparse_partitions[0]
        .seq;

      assert_eq!(
        vec![0, 2, 3, 5, 6],
        seq_info.fitch.variable.keys().copied().collect_vec()
      );

      assert_eq!(&"~C~~C~~TGTATTGAC".chars().collect_vec(), &seq_info.sequence);
    }

    fitch_forward(&graph, &partitions);

    {
      let seq_info = &graph
        .get_exactly_one_root()?
        .read_arc()
        .payload()
        .read_arc()
        .sparse_partitions[0]
        .seq;

      assert_eq!(&aln[0].seq.chars().collect_vec(), &seq_info.sequence);

      let actual_muts: BTreeMap<_, _> = graph
        .get_edges()
        .iter()
        .enumerate()
        .map(|(i, e)| {
          let get_name = |node_id| {
            let node = graph.get_node(node_id).unwrap().read_arc();
            node.payload().read_arc().name().unwrap().as_ref().to_owned()
          };

          let src = get_name(e.read_arc().source());
          let tar = get_name(e.read_arc().target());

          (
            format!("{src}->{tar}"),
            e.read_arc().payload().read_arc().sparse_partitions[0]
              .subs
              .iter()
              .map(ToString::to_string)
              .collect_vec(),
          )
        })
        .collect();

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

      let actual_indels: BTreeMap<_, _> = graph
        .get_edges()
        .iter()
        .map(|e| {
          (
            e.read_arc().key(),
            e.read_arc().payload().read_arc().sparse_partitions[0]
              .indels
              .iter()
              .map(ToString::to_string)
              .collect_vec(),
          )
        })
        .collect();

      let expected_indels = btreemap! {
        0 => vec!["12--13: T -> -", "14--16: -- -> AC"],
        1 => vec![],
        2 => vec!["11--12: T -> -"],
        3 => vec![],
        4 => vec![],
        5 => vec![],
      };

      assert_eq!(
        json_write_str(&expected_indels, JsonPretty(true))?,
        json_write_str(&actual_indels, JsonPretty(true))?
      );
    }

    Ok(())
  }

  #[test]
  fn test_fitch_seq_info_1() -> Result<(), Report> {
    let seq = "GCATCCCTGTA-NG--".chars().collect_vec();
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let actual = SparseSeqNode::new(&seq, &alphabet)?;

    let expected = indoc! { /* language=json */ r#"{
      "seq": {
        "unknown": [
          [
            12,
            13
          ]
        ],
        "gaps": [
          [
            11,
            12
          ],
          [
            14,
            16
          ]
        ],
        "non_char": [
          [
            11,
            13
          ],
          [
            14,
            16
          ]
        ],
        "composition": {
          "counts": {
            "-": 0,
            "A": 0,
            "B": 0,
            "C": 0,
            "D": 0,
            "G": 0,
            "H": 0,
            "K": 0,
            "M": 0,
            "N": 0,
            "R": 0,
            "S": 0,
            "T": 0,
            "V": 0,
            "W": 0,
            "Y": 0
          },
          "gap": "-"
        },
        "sequence": "GCATCCCTGTA-NG--",
        "fitch": {
          "variable": {},
          "variable_indel": {},
          "composition": {
            "counts": {
              "-": 0,
              "A": 0,
              "B": 0,
              "C": 0,
              "D": 0,
              "G": 0,
              "H": 0,
              "K": 0,
              "M": 0,
              "N": 0,
              "R": 0,
              "S": 0,
              "T": 0,
              "V": 0,
              "W": 0,
              "Y": 0
            },
            "gap": "-"
          }
        }
      },
      "profile": {
        "variable": {},
        "variable_indel": {},
        "fixed": {},
        "fixed_counts": {
          "counts": {
            "-": 0,
            "A": 0,
            "B": 0,
            "C": 0,
            "D": 0,
            "G": 0,
            "H": 0,
            "K": 0,
            "M": 0,
            "N": 0,
            "R": 0,
            "S": 0,
            "T": 0,
            "V": 0,
            "W": 0,
            "Y": 0
          },
          "gap": "-"
        },
        "log_lh": 0.0
      }
    }"#};

    assert_eq!(expected, json_write_str(&actual, JsonPretty(true))?);
    Ok(())
  }

  #[test]
  fn test_fitch_seq_info_2() -> Result<(), Report> {
    let seq = "TCGGCCGTGTRTTG--".chars().collect_vec();
    let alphabet = Alphabet::new(AlphabetName::Nuc, false)?;
    let actual = SparseSeqNode::new(&seq, &alphabet)?;

    let expected = indoc! { /* language=json */ r#"{
      "seq": {
        "unknown": [],
        "gaps": [
          [
            14,
            16
          ]
        ],
        "non_char": [
          [
            14,
            16
          ]
        ],
        "composition": {
          "counts": {
            "-": 0,
            "A": 0,
            "B": 0,
            "C": 0,
            "D": 0,
            "G": 0,
            "H": 0,
            "K": 0,
            "M": 0,
            "N": 0,
            "R": 0,
            "S": 0,
            "T": 0,
            "V": 0,
            "W": 0,
            "Y": 0
          },
          "gap": "-"
        },
        "sequence": "TCGGCCGTGTRTTG--",
        "fitch": {
          "variable": {
            "10": [
              "A",
              "G"
            ]
          },
          "variable_indel": {},
          "composition": {
            "counts": {
              "-": 0,
              "A": 0,
              "B": 0,
              "C": 0,
              "D": 0,
              "G": 0,
              "H": 0,
              "K": 0,
              "M": 0,
              "N": 0,
              "R": 0,
              "S": 0,
              "T": 0,
              "V": 0,
              "W": 0,
              "Y": 0
            },
            "gap": "-"
          }
        }
      },
      "profile": {
        "variable": {},
        "variable_indel": {},
        "fixed": {},
        "fixed_counts": {
          "counts": {
            "-": 0,
            "A": 0,
            "B": 0,
            "C": 0,
            "D": 0,
            "G": 0,
            "H": 0,
            "K": 0,
            "M": 0,
            "N": 0,
            "R": 0,
            "S": 0,
            "T": 0,
            "V": 0,
            "W": 0,
            "Y": 0
          },
          "gap": "-"
        },
        "log_lh": 0.0
      }
    }"#};

    assert_eq!(expected, json_write_str(&actual, JsonPretty(true))?);
    Ok(())
  }
}
