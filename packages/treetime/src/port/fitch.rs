#![allow(dead_code)]
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::gtr::gtr::GTR;
use crate::io::fasta::FastaRecord;
use crate::port::composition::Composition;
use crate::port::constants::{FILL_CHAR, GAP_CHAR, NON_CHAR, VARIABLE};
use crate::port::mutation::{InDel, Mut};
use crate::port::seq_partitions::SeqPartition;
use crate::port::seq_sparse::{
  Deletion, SparseGraph, SparseNode, SparseSeqDis, SparseSeqEdge, SparseSeqInfo, SparseSeqNode, VarPos,
};
use crate::seq::range::range_contains;
use crate::seq::range_complement::range_complement;
use crate::seq::range_difference::range_difference;
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

#[derive(Debug)]
pub struct PartitionModel<'g> {
  aln: Vec<FastaRecord>,
  gtr: &'g GTR,
}

fn attach_seqs_to_graph<'g>(graph: &SparseGraph<'g>, partitions: &[PartitionModel<'g>]) -> Result<(), Report> {
  graph.data().write_arc().sparse_partitions = partitions
    .iter()
    .map(|PartitionModel { aln, gtr }| {
      let length = get_common_length(aln)?;
      Ok(SeqPartition::new(gtr, length))
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

        SparseSeqNode::new(&sequence, gtr.alphabet())
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
      let &SeqPartition { gtr, length } = &sparse_partitions[si];

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
              None => gtr.alphabet().get_profile(child.sequence[pos]).view(),
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
          sequence[pos] = gtr.alphabet().char(i);
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
          .filter(|&x| gtr.alphabet().contains(x))
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
            let child_profiles = states
              .iter()
              .map(|&c| gtr.alphabet().get_profile(c).view())
              .collect_vec();
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
          composition: Composition::default(),
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

fn fitch_forward(graph: &mut SparseGraph) {
  let sparse_partitions = &graph.data().read_arc().sparse_partitions;
  let n_partitions = sparse_partitions.len();

  graph.par_iter_breadth_first_forward(|mut node| {
    #[allow(clippy::needless_range_loop)]
    for si in 0..n_partitions {
      let &SeqPartition { gtr, length } = &sparse_partitions[si];
      let SparseSeqInfo {
        gaps,
        sequence,
        composition,
        non_char,
        fitch: SparseSeqDis {
          variable,
          variable_indel,
          ..
        },
        ..
      } = &mut node.payload.sparse_partitions[si].seq;

      if node.is_root {
        for (pos, p) in variable {
          let i = p.dis.argmax().unwrap();
          let state = gtr.alphabet().char(i);
          p.state = Some(state);
          sequence[*pos] = state;
        }
        // process indels as majority rule at the root
        for (r, indel) in variable_indel.iter() {
          if indel.deleted > indel.ins {
            gaps.push(*r);
          }
        }
        for r in gaps.iter() {
          sequence[r.0..r.1].fill(GAP_CHAR);
        }
        *composition = Composition::with_sequence(sequence.iter().copied(), gtr.alphabet().iter().copied());
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
        for (pos, p) in variable.iter_mut() {
          let pnuc = parent.sequence[*pos];
          // check whether parent is in child profile (sum>0 --> parent state is in profile)
          let parent_in_profile = (gtr.alphabet().get_profile(pnuc) * &p.dis).sum() > 0.0;
          if parent_in_profile {
            sequence[*pos] = pnuc;
          } else {
            let i = p.dis.argmax().unwrap();
            let cnuc = gtr.alphabet().char(i);
            sequence[*pos] = cnuc;
            let m = Mut {
              pos: *pos,
              qry: cnuc,
              reff: pnuc,
            };
            composition.add_mutation(&m);
            edge.muts.push(m);
          }
          p.state = Some(sequence[*pos]);
        }

        for (&pos, pvar) in &parent.fitch.variable {
          if variable.contains_key(&pos) {
            continue;
          }

          // NOTE: access to full_seq would not be necessary if we had saved the
          // child state of variable positions in the backward pass
          let node_nuc = sequence[pos];
          if let Some(pvar_state) = pvar.state {
            if pvar_state != node_nuc {
              let m = Mut {
                pos,
                qry: node_nuc,
                reff: pvar_state,
              };
              composition.add_mutation(&m);
              edge.muts.push(m);
            }
          }
        }

        // Process indels. Gaps where the children disagree, need to be decided by also looking at parent.
        for (r, indel) in variable_indel.iter() {
          let gap_in_parent = if parent.gaps.contains(r) { 1 } else { 0 };
          if indel.deleted + gap_in_parent > indel.ins {
            gaps.push(*r);
            if gap_in_parent == 0 {
              // If the gap is not in parent, add deletion.
              let indel = InDel {
                range: *r,
                seq: sequence[r.0..r.1].to_owned(),
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
              deletion: true,
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
            deletion: true,
          };
          composition.add_indel(&indel);
          edge.indels.push(indel);
        }
      }
    }
    GraphTraversalContinuation::Continue
  });
}

fn fitch_cleanup() {}

pub fn compress_sequences<'g>(graph: &mut SparseGraph<'g>, partitions: &'g [PartitionModel<'g>]) -> Result<(), Report> {
  attach_seqs_to_graph(graph, partitions)?;
  fitch_backwards(graph);
  fitch_forward(graph);
  fitch_cleanup();
  Ok(())
}

/// Reconstruct ancestral sequences using Fitch parsimony.
///
/// Calls visitor function for every ancestral node, providing the node itself and its reconstructed sequence.
/// Optionally reconstructs leaf sequences.
pub fn ancestral_reconstruction_fitch(
  graph: &SparseGraph,
  include_leaves: bool,
  mut visitor: impl FnMut(&SparseNode, Vec<char>),
) -> Result<(), Report> {
  let sparse_partitions = &graph.data().read_arc().sparse_partitions;
  let n_partitions = sparse_partitions.len();

  graph.iter_depth_first_preorder_forward(|node| {
    if node.is_root {
      return;
    }
    if !include_leaves && node.is_leaf {
      return;
    }

    let seq = (0..n_partitions)
      .flat_map(|si| {
        let (parent, edge) = node.get_exactly_one_parent().unwrap();
        let parent = &parent.read_arc().sparse_partitions[si];
        let edge = &edge.read_arc().sparse_partitions[si];

        let node = &node.payload.sparse_partitions[si].seq;

        let mut seq = parent.seq.sequence.clone();

        // Implant mutations
        for m in &edge.muts {
          seq[m.pos] = m.qry;
        }

        // Implant indels
        for indel in &edge.indels {
          if indel.deletion {
            seq[indel.range.0..indel.range.1].fill(GAP_CHAR);
          } else {
            seq[indel.range.0..indel.range.1].copy_from_slice(&indel.seq);
          }
        }

        // At the node itself, mask whatever is unknown in the node.
        for r in &node.unknown {
          let ambig_char = sparse_partitions[si].gtr.alphabet.ambiguous();
          seq[r.0..r.1].fill(ambig_char);
        }

        for (pos, p) in &node.fitch.variable {
          seq[*pos] = sparse_partitions[si].code(&p.dis);
        }

        seq
      })
      .collect();

    visitor(&node.payload, seq);
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
  use crate::gtr::get_gtr::{jc69, JC69Params};
  use crate::io::fasta::read_many_fasta_str;
  use crate::io::json::{json_write_str, JsonPretty};
  use crate::io::nwk::nwk_read_str;
  use crate::utils::string::vec_to_string;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;

  #[test]
  fn test_seq_ancestral_reconstruction() -> Result<(), Report> {
    rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;

    let inputs = read_many_fasta_str(indoc! {r#"
      >A
      ACATCGCCNNA--G
      >B
      GCATCCCTGTA-NG
      >C
      CCGGCGATGTATTG
      >D
      TCGGCCGTGTRTTG
    "#})?;

    let expected = read_many_fasta_str(indoc! {r#"
      >AB
      GCATCGCTGTATTG
      >CD
      TCGGCGGTGTATTG
      >root
      TCGGCGGTGTATTG
    "#})?;

    let mut graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let gtr = &jc69(JC69Params::default())?;
    let partitions = vec![PartitionModel { gtr, aln: inputs }];
    compress_sequences(&mut graph, &partitions).unwrap();

    let mut actual = BTreeMap::new();
    ancestral_reconstruction_fitch(&graph, false, |node, seq| {
      actual.insert(node.name.clone(), vec_to_string(seq));
    })?;

    assert_eq!(
      json_write_str(&expected, JsonPretty(false))?,
      json_write_str(&actual, JsonPretty(false))?
    );

    Ok(())
  }

  // #[test]
  // fn test_seq_ancestral_reconstruction2() -> Result<(), Report> {
  //   rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
  //
  //   let inputs = read_many_fasta_str(indoc! {r#"
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
  //   "#})?;
  //
  //   let expected = read_many_fasta_str(indoc! {r#"
  //     >root
  //     ACAGCCATGTATTG--
  //     >AB
  //     ACATCCCTGTA-TG--
  //     >CD
  //     CCGGCCATGTATTG--
  //   "#})?;
  //
  //   let mut graph: SparseGraph = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
  //
  //   let gtr = &jc69(JC69Params::default())?;
  //   let partitions = vec![PartitionModel { gtr, aln: inputs }];
  //   compress_sequences(&mut graph, &partitions).unwrap();
  //
  //   let mut actual = BTreeMap::new();
  //   ancestral_reconstruction_fitch(&graph, false, |node, seq| {
  //     actual.insert(node.name.clone(), vec_to_string(seq));
  //   })?;
  //
  //   assert_eq!(
  //     json_write_str(&expected, JsonPretty(false))?,
  //     json_write_str(&actual, JsonPretty(false))?
  //   );
  //
  //   Ok(())
  // }

  // #[test]
  // fn test_seq_ancestral_reconstruction_with_leaves() -> Result<(), Report> {
  //   rayon::ThreadPoolBuilder::new().num_threads(1).build_global()?;
  //
  //   let inputs = BTreeMap::from([
  //     (o!("A"), o!("ACATCGCCNNA--G")),
  //     (o!("B"), o!("GCATCCCTGTA-NG")),
  //     (o!("C"), o!("CCGGCGATGTATTG")),
  //     (o!("D"), o!("TCGGCCGTGTRTTG")),
  //   ]);
  //
  //   #[rustfmt::skip]
  //   let expected = BTreeMap::from([
  //     (o!("A"),    o!("ACATCGCCNNA--G")),
  //     (o!("AB"),   o!("GCATCGCTGTATTG")),
  //     (o!("B"),    o!("GCATCCCTGTA-NG")),
  //     (o!("C"),    o!("CCGGCGATGTATTG")),
  //     (o!("CD"),   o!("TCGGCGGTGTATTG")),
  //     (o!("D"),    o!("TCGGCCGTGTRTTG")),
  //     (o!("root"), o!("TCGGCGGTGTATTG")),
  //   ]);
  //
  //   let graph = nwk_read_str::<Node, Edge, ()>("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
  //
  //   compress_sequences(&inputs, &graph).unwrap();
  //
  //   let mut actual = BTreeMap::new();
  //   ancestral_reconstruction_fitch(&graph, true, |node, seq| {
  //     actual.insert(node.name.clone(), vec_to_string(seq.to_owned()));
  //   })?;
  //
  //   assert_eq!(
  //     json_write_str(&expected, JsonPretty(false))?,
  //     json_write_str(&actual, JsonPretty(false))?
  //   );
  //
  //   Ok(())
  // }
}
