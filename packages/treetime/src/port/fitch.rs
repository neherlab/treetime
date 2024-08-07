#![allow(dead_code)]
use crate::alphabet::alphabet::Alphabet;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::gtr::gtr::GTR;
use crate::io::fasta::FastaRecord;
use crate::port::seq_partitions::SeqPartition;
use crate::port::seq_sparse::{
  Deletion, SparseGraph, SparseNode, SparseSeqDis, SparseSeqEdge, SparseSeqInfo, SparseSeqNode, VarPos,
};
use crate::seq::range::range_contains;
use crate::seq::range_complement::range_complement;
use crate::seq::range_intersection::{range_intersection, range_intersection_iter};
use crate::utils::manyzip::Manyzip;
use crate::utils::ndarray::product_axis;
use crate::{make_error, make_internal_report, make_report};
use approx::UlpsEq;
use eyre::{Report, WrapErr};
use itertools::{izip, Itertools};
use maplit::btreemap;
use ndarray::{stack, Array1, Array2, Axis};
use ndarray_stats::QuantileExt;

const NON_CHAR: char = '.';
const VARIABLE: char = '~';
const FILL_CHAR: char = ' ';
const GAP_CHAR: char = '-';

#[derive(Debug)]
pub struct PartitionModel<'g> {
  aln: Vec<FastaRecord>,
  gtr: &'g GTR,
}

fn attach_seqs_to_graph<'a, 'g>(
  graph: &SparseGraph<'a, 'g>,
  alphabet: &'a Alphabet,
  partitions: &[PartitionModel<'g>],
) -> Result<(), Report> {
  graph.data().write_arc().sparse_partitions = partitions
    .iter()
    .map(|PartitionModel { aln, gtr }| {
      let length = get_common_length(aln)?;
      Ok(SeqPartition::new(gtr, length, alphabet))
    })
    .collect::<Result<_, Report>>()?;

  for leaf in graph.get_leaves() {
    let mut leaf = leaf.read_arc().payload().write_arc();

    let leaf_name = leaf.name.as_ref().ok_or_else(|| {
      make_report!("Expected all leaf nodes to have names, such that they can be matched to their corresponding sequences. But found a leaf node that has no name.")
    })?.to_owned();

    let sparse_partitions = partitions
      .iter()
      .map(|PartitionModel { aln, gtr }| {
        // TODO(perf): this might be slow if there are many sequences
        let leaf_fasta = aln
          .iter()
          .find(|fasta| fasta.seq_name == leaf_name)
          .ok_or_else(|| make_internal_report!("Leaf sequence not found: '{leaf_name}'"))?;

        // TODO(perf): unnecessary copy of sequence data. Neither String, nor &[char] works well for us, it seems.
        // We probably want a custom class for sequences. Sequences should be instantiated in the fasta parser and
        // never need a copy like here.
        let sequence = leaf_fasta.seq.chars().collect::<Vec<_>>();

        SparseSeqNode::new(&sequence, alphabet)
      })
      .collect::<Result<_, Report>>()?;

    *leaf = SparseNode {
      name: Some(leaf_name),
      sparse_partitions,
    }
  }

  Ok(())
}

fn fitch_backwards(graph: &mut SparseGraph) {
  let sparse_partitions = &graph.data().read_arc().sparse_partitions;
  let n_partitions = sparse_partitions.len();

  graph.par_iter_breadth_first_backward(|mut node| {
    if graph.is_leaf(node.key) {
      return GraphTraversalContinuation::Continue;
    }

    for (_, edge) in &node.children {
      edge.write_arc().sparse_partitions = vec![SparseSeqEdge::default(); n_partitions];
    }

    #[allow(clippy::needless_range_loop)]
    for si in 0..n_partitions {
      let &SeqPartition { gtr, length, alphabet } = &sparse_partitions[si];

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

      let mut gaps = range_intersection_iter(children.iter().map(|(c, _)| &c.gaps)).collect_vec();
      let unknown = range_intersection_iter(children.iter().map(|(c, _)| &c.unknown)).collect_vec();
      let non_char = range_intersection_iter(children.iter().map(|(c, _)| &c.non_char)).collect_vec();
      let non_gap = range_complement(&[(0, length)], &[gaps.clone()]); // FIXME(perf): unnecessary clone

      let mut seq_dis = SparseSeqDis::default();
      let mut sequence = vec![FILL_CHAR; length];

      for r in &non_char {
        sequence[r.0..r.1].fill(NON_CHAR);
      }

      // Process all positions where the children are variable.
      // Need to account for parts of the sequence transmitted along edges.
      let variable_positions = children
        .iter()
        .flat_map(|(c, e)| c.fitch.variable.keys().copied())
        .collect_vec();
      for pos in variable_positions {
        // Collect child profiles (1D vectors)
        let child_profiles = children
          .iter()
          .filter_map(|(child, edge)| {
            if range_contains(&edge.transmission, pos) {
              return None; // transmission field is not currently used
            }
            if range_contains(&child.non_char, pos) {
              return None; // this position does not have character state information
            }
            let state = match child.fitch.variable.get(&pos) {
              Some(var_pos) => var_pos.dis.view(),
              None => alphabet.get_profile(child.sequence[pos]).view(),
            };
            Some(state)
          })
          .collect_vec();

        // Stack child profiles on top of each other into a 2D matrix
        let child_profiles: Array2<f64> = stack(Axis(0), &child_profiles).unwrap();

        // Calculate Fitch parsimony.
        // If we save the states of the children for each position that is variable in the node,
        // then we would not need the full sequences in the forward pass.
        let intersection = product_axis(&child_profiles, Axis(0));
        if intersection.sum().ulps_eq(&1.0, f64::EPSILON, 5) {
          let i = intersection.argmax().unwrap();
          sequence[pos] = alphabet.char(i);
        } else if intersection.sum() > 1.0 {
          seq_dis.variable.insert(pos, VarPos::new(intersection, None));
          sequence[pos] = VARIABLE;
        } else {
          let union = child_profiles.sum_axis(Axis(0));
          let dis = Array1::from_iter(union.iter().map(|&x| if x > 0.0 { 1.0 } else { 0.0 }));
          seq_dis.variable.insert(pos, VarPos::new(dis, None));
          sequence[pos] = VARIABLE;
        }
      }

      // Process all positions where the children are fixed or completely unknown in some children.

      // Gather state sets for each position across child sequences
      // TODO(perf): avoid copying and allocations
      let child_state_sets = Manyzip(children.iter().map(|(c, e)| c.sequence.iter().copied()).collect_vec());

      // Zip these states with node sequence
      let state_zip = izip!(sequence.iter_mut(), child_state_sets.into_iter());
      for (pos, (nuc, child_states)) in state_zip.enumerate() {
        if *nuc != FILL_CHAR {
          continue;
        }

        let determined_states = child_states
          .into_iter()
          .filter(|&x| alphabet.contains(x))
          .unique()
          .collect_vec();

        // Find the state of the current node at this position
        *nuc = match determined_states.as_slice() {
          [state] => {
            // All children have the same state, that will be the state of the current node
            *state
          }
          [] => {
            // No child states. Impossible
            unreachable!("No child states. This is impossible");
          }
          states => {
            // Child states differ. This is variable state.
            // Save child states and postpone the decision until forward pass.
            let child_profiles = states.iter().map(|&c| alphabet.get_profile(c).view()).collect_vec();
            let child_profiles: Array2<f64> = stack(Axis(0), &child_profiles).unwrap();
            let summed_profile = child_profiles.sum_axis(Axis(0));
            seq_dis.variable.insert(pos, VarPos::new(summed_profile, None));
            VARIABLE
          }
        };
      }

      // Process insertions and deletions.

      // 1) seq_info.gaps is the intersection of child gaps, i.e. this is gap if and only if all children have a gap
      //    --> hence we find positions where children differ in terms of gap presence absence by intersecting
      //        the child gaps with the complement of the parent gaps
      for (child, _) in &children {
        for r in range_intersection(&[non_gap.clone(), child.gaps.clone()]) {
          let indel = seq_dis.variable_indel.entry(r).or_insert_with(|| Deletion {
            deleted: 0,
            ins: n_children,
            alt: sequence[r.0..r.1].to_owned(),
          });
          indel.deleted += 1;
          indel.ins -= 1;
        }
      }

      // 2) if a gap is variable in a child and the parent, we need to pull this down to the parent
      for (child, _) in &children {
        for r in child.fitch.variable_indel.keys() {
          if let Some(indel) = seq_dis.variable_indel.get_mut(r) {
            indel.deleted += 1;
            indel.ins -= 1;
          }
        }
      }

      // 3) if all children are compatible with a gap, we add the gap back to the gap collection and remove the
      // variable site (nothing needs doing in the case where all children are compatible with non-gap)
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
          composition: btreemap! {},
        },
        profile: SparseSeqDis::default(),
        msg_to_parents: SparseSeqDis::default(),
        msgs_to_children: btreemap! {},
        msgs_from_children: btreemap! {},
      });
    }

    GraphTraversalContinuation::Continue
  });
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
