#[cfg(test)]
mod tests {
  use super::helpers::{ambiguity, blocks, del, ins, interleaved, layout_result, piece, seq, seqs, sub, subs, whole};
  use crate::partition::assembly::AlignmentLayout;
  use crate::partition::resolve::PartitionLayout;
  use crate::seq::indel::{InDel, compose_indels};
  use crate::seq::mutation::{Mutation, MutationTrack, Sub, combine_edge_mutations};
  use crate::seq::sink::AmbiguousSites;
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use std::slice::from_ref;
  use treetime_primitives::Seq;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::interleaved(interleaved(6),   &["AC", "GT", "CA"], "AGCCTA")]
  #[case::blocks(     blocks(&[2, 4]),  &["AC", "GTCA"],     "ACGTCA")]
  #[case::whole(      whole(6),         &["ACGTCA"],         "ACGTCA")]
  #[trace]
  fn test_assembly_sequence_writes_each_partition_to_its_columns(
    #[case] layout: AlignmentLayout,
    #[case] parts: &[&str],
    #[case] expected: &str,
  ) -> Result<(), Report> {
    assert_eq!(seq(expected), layout.assemble_sequence(&seqs(parts))?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::interleaved(interleaved(6),  vec![vec![(0, 'R')], vec![(0, 'Y')], vec![(1, 'N')]], &[(0, 'R'), (1, 'Y'), (5, 'N')])]
  #[case::blocks(     blocks(&[3, 3]), vec![vec![(0, 'R')], vec![(0, 'Y'), (2, 'K')]],       &[(0, 'R'), (3, 'Y'), (5, 'K')])]
  #[case::whole(      whole(6),        vec![vec![(0, 'R'), (4, 'Y')]],                       &[(0, 'R'), (4, 'Y')]          )]
  #[trace]
  fn test_assembly_ambiguity_maps_same_local_position_to_own_sites(
    #[case] layout: AlignmentLayout,
    #[case] parts: Vec<Vec<(usize, char)>>,
    #[case] expected: &[(usize, char)],
  ) -> Result<(), Report> {
    let parts = parts.iter().map(|sites| ambiguity(sites)).collect::<Vec<AmbiguousSites>>();
    assert_eq!(ambiguity(expected), layout.assemble_ambiguity(&parts)?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::interleaved(interleaved(6),  vec![vec!["C2T"], vec!["A1G"], vec![]], &["A2G", "C4T"])]
  #[case::blocks(     blocks(&[3, 3]), vec![vec!["A3G"], vec!["C1T"]],         &["A3G", "C4T"])]
  #[case::whole(      whole(6),        vec![vec!["A3G", "C4T"]],               &["A3G", "C4T"])]
  #[trace]
  fn test_assembly_substitutions_map_positions_and_sort(
    #[case] layout: AlignmentLayout,
    #[case] parts: Vec<Vec<&str>>,
    #[case] expected: &[&str],
  ) -> Result<(), Report> {
    let parts = parts.iter().map(|part| subs(part)).collect::<Vec<Vec<Sub>>>();
    assert_eq!(subs(expected), layout.assemble_substitutions(&parts)?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::deletion_across_block_border( blocks(&[3, 3]), vec![vec![del(1, 3, "CG")], vec![del(0, 2, "TA")]],                   vec![del(1, 5, "CGTA")]                 )]
  #[case::deletion_across_interleaved(  interleaved(6),  vec![vec![del(1, 2, "T")], vec![del(1, 2, "A")], vec![del(0, 1, "G")]], vec![del(2, 5, "GTA")]                  )]
  #[case::insertion_across_interleaved( interleaved(6),  vec![vec![ins(0, 1, "A")], vec![ins(0, 1, "C")], vec![]],               vec![ins(0, 2, "AC")]                   )]
  #[case::local_run_apart_in_alignment( interleaved(9),  vec![vec![del(0, 2, "AC")], vec![], vec![]],                            vec![del(0, 1, "A"), del(3, 4, "C")]    )]
  #[case::kinds_stay_apart(             blocks(&[2, 2]), vec![vec![del(1, 2, "A")], vec![ins(0, 1, "C")]],                       vec![del(1, 2, "A"), ins(2, 3, "C")]    )]
  #[case::whole_keeps_adjacent_events(  whole(4),        vec![vec![ins(0, 1, "A"), ins(1, 2, "C")]],                             vec![ins(0, 1, "A"), ins(1, 2, "C")]    )]
  #[trace]
  fn test_assembly_indels_join_pass_pieces_into_maximal_runs(
    #[case] layout: AlignmentLayout,
    #[case] parts: Vec<Vec<InDel>>,
    #[case] expected: Vec<InDel>,
  ) -> Result<(), Report> {
    assert_eq!(expected, layout.assemble_indels(&parts)?);
    Ok(())
  }

  #[test]
  fn test_assembly_edge_mutations_put_deletions_before_insertions() -> Result<(), Report> {
    let layout = interleaved(6);
    let substitutions = vec![vec![], subs(&["A1G"]), vec![]];
    let indels = vec![vec![ins(0, 1, "A"), del(1, 2, "C")], vec![del(1, 2, "T")], vec![]];
    let track = MutationTrack::Nucleotide;
    let expected = vec![
      Mutation::substitution(track.clone(), sub("A2G")),
      Mutation::indel(track.clone(), &del(3, 5, "CT"))?,
      Mutation::indel(track.clone(), &ins(0, 1, "A"))?,
    ];
    assert_eq!(
      expected,
      layout.assemble_edge_mutations(&substitutions, &indels, &track)?
    );
    Ok(())
  }

  #[test]
  fn test_assembly_adjacent_composed_insertions_stay_two_events() -> Result<(), Report> {
    let layout = interleaved(6);
    let indels = vec![ins(0, 1, "A"), ins(1, 2, "C")];
    let pieces = layout.split_indels(&indels)?;
    let expected_pieces = vec![
      vec![piece(ins(0, 1, "A"), 0, (0, 1))],
      vec![piece(ins(0, 1, "C"), 1, (1, 2))],
      vec![],
    ];
    assert_eq!(expected_pieces, pieces);
    assert_eq!(indels, layout.assemble_indels(&pieces)?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::insertion_before_deletion( vec![ins(0, 1, "A")], vec![del(3, 5, "CT")])]
  #[case::same_range_parent_first(   vec![del(2, 3, "C")], vec![ins(2, 3, "A")] )]
  #[trace]
  fn test_assembly_composed_edge_mutations_keep_unsplit_order(
    #[case] parent: Vec<InDel>,
    #[case] child: Vec<InDel>,
    #[values(interleaved(6), blocks(&[3, 3]), whole(6))] layout: AlignmentLayout,
  ) -> Result<(), Report> {
    let track = MutationTrack::Nucleotide;
    let expected = combine_edge_mutations(vec![], &compose_indels(&parent, &child), &track)?;
    let pieces = layout.compose_indels(&layout.split_indels(&parent)?, &layout.split_indels(&child)?)?;
    let substitutions = vec![vec![]; pieces.len()];
    assert_eq!(expected, layout.assemble_edge_mutations(&substitutions, &pieces, &track)?);
    Ok(())
  }

  #[test]
  fn test_assembly_rejects_pass_pieces_next_to_split_pieces() {
    let layout = interleaved(6);
    let pieces = vec![vec![del(0, 1, "A")], vec![piece(ins(0, 1, "C"), 0, (1, 2))], vec![]];
    assert_error!(
      layout.assemble_indels(&pieces),
      "An edge of alignment 1 holds indel pieces of the marginal passes next to pieces split from composed indels. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_assembly_compose_and_split_back_across_borders() -> Result<(), Report> {
    let layout = blocks(&[3, 3, 2]);
    let parent = vec![vec![del(2, 3, "G")], vec![del(0, 1, "T")], vec![]];
    let child = vec![vec![], vec![ins(2, 3, "C")], vec![ins(0, 1, "A")]];
    let expected = vec![
      vec![piece(del(2, 3, "G"), 0, (2, 4))],
      vec![piece(del(0, 1, "T"), 0, (2, 4)), piece(ins(2, 3, "C"), 1, (5, 7))],
      vec![piece(ins(0, 1, "A"), 1, (5, 7))],
    ];
    let composed = layout.compose_indels(&parent, &child)?;
    assert_eq!(expected, composed);
    assert_eq!(
      vec![del(2, 4, "GT"), ins(5, 7, "CA")],
      layout.assemble_indels(&composed)?
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::adjacent_insertions(        vec![ins(0, 1, "A")],  vec![ins(1, 2, "C")],  vec![ins(0, 1, "A"), ins(1, 2, "C")]  )]
  #[case::deletion_and_insertion(     vec![del(2, 4, "GT")], vec![ins(5, 7, "CA")], vec![del(2, 4, "GT"), ins(5, 7, "CA")] )]
  #[case::adjacent_deletions_merge(   vec![del(1, 3, "CG")], vec![del(3, 5, "TA")], vec![del(1, 5, "CGTA")]                )]
  #[case::insertion_inside_deletion(  vec![del(1, 6, "CGTAC")], vec![ins(3, 4, "G")], vec![del(1, 6, "CGTAC"), ins(3, 4, "G")] )]
  #[trace]
  fn test_assembly_split_of_composed_indels_assembles_back(
    #[case] parent: Vec<InDel>,
    #[case] child: Vec<InDel>,
    #[case] expected: Vec<InDel>,
    #[values(interleaved(9), blocks(&[3, 3, 3]), blocks(&[1, 1, 7]), whole(9))] layout: AlignmentLayout,
  ) -> Result<(), Report> {
    let composed = compose_indels(&parent, &child);
    assert_eq!(expected, composed);
    assert_eq!(composed, layout.assemble_indels(&layout.split_indels(&composed)?)?);
    let composed_by_layout = layout.compose_indels(&layout.split_indels(&parent)?, &layout.split_indels(&child)?)?;
    assert_eq!(composed, layout.assemble_indels(&composed_by_layout)?);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::no_partition( 3, vec![],                       "An alignment layout needs at least one partition. This is an internal error. Please report it to developers.")]
  #[case::shared(       3, vec![vec![0, 1], vec![1, 2]],    "partition 1 and partition 2 both take column 1 of alignment 1. This is an internal error. Please report it to developers.")]
  #[case::uncovered(    3, vec![vec![0], vec![2]],          "No partition takes column 1 of alignment 1. This is an internal error. Please report it to developers.")]
  #[case::outside(      2, vec![vec![0, 1, 2]],             "partition 1 takes column 2, but alignment 1 has 2 columns. This is an internal error. Please report it to developers.")]
  #[case::descending(   3, vec![vec![2, 0], vec![1]],       "Columns of partition 1 are not ascending: 2 before 0. This is an internal error. Please report it to developers.")]
  #[trace]
  fn test_assembly_layout_rejects_columns_that_do_not_cover_alignment_once(
    #[case] length: usize,
    #[case] partitions: Vec<Vec<usize>>,
    #[case] message: &str,
  ) {
    let partitions = partitions.iter().map(Vec::as_slice).collect::<Vec<_>>();
    assert_error!(layout_result(length, &partitions), message);
  }

  #[test]
  fn test_assembly_layout_rejects_partitions_of_two_alignments() {
    let first = PartitionLayout {
      number: 1,
      name: None,
      alignment: 1,
      columns: vec![0],
    };
    let second = PartitionLayout {
      number: 2,
      name: Some("cp2".to_owned()),
      alignment: 2,
      columns: vec![1],
    };
    assert_error!(
      AlignmentLayout::new(2, [&first, &second]),
      "partition 2 'cp2' belongs to alignment 2, but the layout is of alignment 1. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_assembly_rejects_pieces_that_leave_out_part_of_their_indel() {
    let layout = interleaved(6);
    let pieces = vec![
      vec![piece(del(1, 2, "T"), 0, (2, 5))],
      vec![],
      vec![piece(del(0, 1, "G"), 0, (2, 5))],
    ];
    assert_error!(
      layout.assemble_indels(&pieces),
      "The pieces of indel 2--5 (Deletion) of alignment 1 cover the positions 2, 3. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_assembly_rejects_piece_outside_its_partition() {
    let layout = interleaved(6);
    let pieces = vec![vec![del(1, 3, "TA")], vec![], vec![]];
    assert_error!(
      layout.assemble_indels(&pieces),
      "Indel 1--3: TA -> -- does not fit partition 1, which has 2 columns. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_assembly_rejects_wrong_number_of_partition_values() {
    let layout = interleaved(6);
    assert_error!(
      layout.assemble_sequence(&seqs(&["AC", "GT"])),
      "Alignment 1 has 3 partitions, but assembly received 2 sequences. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_assembly_split_rejects_piece_of_split_indel() {
    let layout = interleaved(6);
    assert_error!(
      layout.split_indels(&[piece(del(1, 2, "T"), 0, (2, 5))]),
      "Indel 1--2: T -> - of alignment 1 is a piece of the indel at 2--5. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_assembly_split_rejects_indel_outside_alignment() {
    let layout = interleaved(6);
    assert_error!(
      layout.split_indels(&[del(5, 7, "TA")]),
      "Indel 5--7: TA -> -- does not fit alignment 1, which has 6 columns. This is an internal error. Please report it to developers."
    );
  }

  proptest! {
    #[test]
    fn test_prop_assembly_split_assembles_back((layout, indels) in generators::layout_and_indels()) {
      let pieces = layout.split_indels(&indels).unwrap();
      prop_assert_eq!(indels, layout.assemble_indels(&pieces).unwrap());
    }

    #[test]
    fn test_prop_assembly_whole_layout_passes_values_through(
      (length, indels) in generators::length_and_indels(),
      chars in "[ACGTRYN]{24}",
    ) {
      let layout = whole(length);
      let pieces = layout.split_indels(&indels).unwrap();
      prop_assert_eq!(&vec![indels.clone()], &pieces);
      prop_assert_eq!(&indels, &layout.assemble_indels(&pieces).unwrap());
      let sequence = Seq::from(&seq(&chars).as_slice()[..length]);
      prop_assert_eq!(&sequence, &layout.assemble_sequence(from_ref(&sequence)).unwrap());
    }

    #[test]
    fn test_prop_assembly_pass_pieces_give_unsplit_indels(
      (layout, columns, classes, chars) in generators::layout_and_gap_classes(),
    ) {
      let expected = generators::maximal_runs(&classes, &chars);
      let pieces = columns
        .iter()
        .map(|columns| {
          let local_classes = columns.iter().map(|&column| classes[column]).collect::<Vec<_>>();
          let local_chars: Seq = columns.iter().map(|&column| chars[column]).collect();
          generators::maximal_runs(&local_classes, &local_chars)
        })
        .collect::<Vec<_>>();
      prop_assert_eq!(expected, layout.assemble_indels(&pieces).unwrap());
    }
  }

  mod generators {
    use super::super::helpers::layout;
    use crate::partition::assembly::AlignmentLayout;
    use crate::seq::indel::{InDel, InDelKind};
    use itertools::Itertools;
    use proptest::collection::vec;
    use proptest::prelude::*;
    use treetime_primitives::Seq;

    pub(super) fn layout_and_indels() -> impl Strategy<Value = (AlignmentLayout, Vec<InDel>)> {
      columns_by_partition()
        .prop_flat_map(|(length, columns)| (Just(columns), indels(length)))
        .prop_map(|(columns, indels)| {
          let length = columns.iter().map(Vec::len).sum();
          (layout(length, &columns.iter().map(Vec::as_slice).collect_vec()), indels)
        })
    }

    pub(super) fn length_and_indels() -> impl Strategy<Value = (usize, Vec<InDel>)> {
      (1_usize..=24).prop_flat_map(|length| (Just(length), indels(length)))
    }

    pub(super) fn layout_and_gap_classes()
    -> impl Strategy<Value = (AlignmentLayout, Vec<Vec<usize>>, Vec<Option<InDelKind>>, Seq)> {
      columns_by_partition().prop_flat_map(|(length, columns)| {
        let classes = vec(
          prop_oneof![
            Just(None),
            Just(Some(InDelKind::Deletion)),
            Just(Some(InDelKind::Insertion))
          ],
          length,
        );
        let chars = proptest::string::string_regex(&format!("[ACGT]{{{length}}}")).unwrap();
        (Just(columns), classes, chars).prop_map(move |(columns, classes, chars)| {
          let layout = layout(length, &columns.iter().map(Vec::as_slice).collect_vec());
          (layout, columns, classes, Seq::try_from_str(&chars).unwrap())
        })
      })
    }

    pub(super) fn maximal_runs(classes: &[Option<InDelKind>], chars: &Seq) -> Vec<InDel> {
      let mut runs: Vec<InDel> = vec![];
      for (position, (class, &state)) in classes.iter().zip(chars.iter()).enumerate() {
        let Some(kind) = *class else { continue };
        match runs
          .iter_mut()
          .rev()
          .find(|run| run.kind == kind && run.range.1 == position)
        {
          Some(run) => {
            run.range.1 += 1;
            run.seq.push(state);
          },
          None => runs.push(InDel {
            range: (position, position + 1),
            seq: Seq::from(&[state][..]),
            kind,
            split_from: None,
          }),
        }
      }
      runs.sort_by_key(|run| (run.kind == InDelKind::Insertion, run.range));
      runs
    }

    fn columns_by_partition() -> impl Strategy<Value = (usize, Vec<Vec<usize>>)> {
      (1_usize..=24)
        .prop_flat_map(|length| vec(0_usize..4, length))
        .prop_map(|owners| {
          let columns = (0..4)
            .map(|owner| owners.iter().positions(|&o| o == owner).collect_vec())
            .filter(|columns| !columns.is_empty())
            .collect_vec();
          (owners.len(), columns)
        })
    }

    fn indels(length: usize) -> impl Strategy<Value = Vec<InDel>> {
      vec(
        (0..length, any::<bool>(), "[ACGT]{1,5}").prop_map(move |(lo, deletion, chars)| {
          let hi = (lo + chars.len()).min(length);
          InDel {
            range: (lo, hi),
            seq: Seq::from(&Seq::try_from_str(&chars).unwrap().as_slice()[..hi - lo]),
            kind: if deletion {
              InDelKind::Deletion
            } else {
              InDelKind::Insertion
            },
            split_from: None,
          }
        }),
        0..8,
      )
    }
  }
}

mod helpers {
  use crate::partition::assembly::AlignmentLayout;
  use crate::partition::resolve::PartitionLayout;
  use crate::seq::indel::{InDel, SplitFrom};
  use crate::seq::mutation::Sub;
  use crate::seq::sink::AmbiguousSites;
  use crate::test_utils::{deletion, insertion};
  use eyre::Report;
  use itertools::{Itertools, izip};
  use treetime_primitives::{AsciiChar, Seq};

  pub(super) fn layout_result(length: usize, partitions: &[&[usize]]) -> Result<AlignmentLayout, Report> {
    let layouts = izip!(1.., partitions)
      .map(|(number, columns)| PartitionLayout {
        number,
        name: None,
        alignment: 1,
        columns: columns.to_vec(),
      })
      .collect_vec();
    AlignmentLayout::new(length, &layouts)
  }

  pub(super) fn layout(length: usize, partitions: &[&[usize]]) -> AlignmentLayout {
    layout_result(length, partitions).unwrap()
  }

  pub(super) fn interleaved(length: usize) -> AlignmentLayout {
    let columns = (0..3)
      .map(|offset| (offset..length).step_by(3).collect_vec())
      .collect_vec();
    layout(length, &columns.iter().map(Vec::as_slice).collect_vec())
  }

  pub(super) fn blocks(sizes: &[usize]) -> AlignmentLayout {
    let columns = sizes
      .iter()
      .scan(0, |start, &size| {
        let block = (*start..*start + size).collect_vec();
        *start += size;
        Some(block)
      })
      .collect_vec();
    layout(sizes.iter().sum(), &columns.iter().map(Vec::as_slice).collect_vec())
  }

  pub(super) fn whole(length: usize) -> AlignmentLayout {
    layout(length, &[&(0..length).collect_vec()])
  }

  pub(super) fn seq(text: &str) -> Seq {
    Seq::try_from_str(text).unwrap()
  }

  pub(super) fn seqs(texts: &[&str]) -> Vec<Seq> {
    texts.iter().map(|text| seq(text)).collect()
  }

  pub(super) fn ambiguity(sites: &[(usize, char)]) -> AmbiguousSites {
    sites
      .iter()
      .map(|&(position, code)| (position, AsciiChar::try_from_char(code).unwrap()))
      .collect()
  }

  pub(super) fn sub(text: &str) -> Sub {
    text.parse().unwrap()
  }

  pub(super) fn subs(texts: &[&str]) -> Vec<Sub> {
    texts.iter().map(|text| sub(text)).collect()
  }

  pub(super) fn del(lo: usize, hi: usize, text: &str) -> InDel {
    deletion((lo, hi), seq(text))
  }

  pub(super) fn ins(lo: usize, hi: usize, text: &str) -> InDel {
    insertion((lo, hi), seq(text))
  }

  pub(super) fn piece(indel: InDel, index: usize, range: (usize, usize)) -> InDel {
    InDel {
      split_from: Some(SplitFrom { index, range }),
      ..indel
    }
  }
}
