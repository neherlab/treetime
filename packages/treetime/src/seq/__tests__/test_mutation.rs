#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::Alphabet;
  use crate::seq::mutation::{
    AlignedMutation, Mutation, MutationEvent, MutationTrack, PositionBase, Sub, compose_substitutions,
    edge_state_change_counts, mutation_event_strings,
  };
  use eyre::Report;
  use treetime_graph::keymap;

  use pretty_assertions::assert_eq;
  use proptest::prelude::*;
  use rstest::rstest;
  use treetime_graph::graph::Graph;
  use treetime_primitives::AsciiChar;
  use treetime_primitives::seq;
  use treetime_utils::assert_error;

  #[test]
  fn test_mutation_sub_display_uses_one_based_position() -> Result<(), Report> {
    let sub = Sub::new(AsciiChar::try_new(b'A')?, 9_usize, AsciiChar::try_new(b'G')?)?;
    assert_eq!("A10G", sub.to_string());
    Ok(())
  }

  #[test]
  fn test_mutation_sub_display_roundtrips_through_from_str() -> Result<(), Report> {
    let expected: Sub = "C123T".parse()?;
    let actual: Sub = expected.to_string().parse()?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_both_empty() -> Result<(), Report> {
    let result = compose_substitutions(&[], &[])?;
    assert_eq!(result, vec![]);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_parent_empty() -> Result<(), Report> {
    let child = vec![helpers::sub(b'A', 0, b'T'), helpers::sub(b'G', 5, b'C')];
    let result = compose_substitutions(&[], &child)?;
    assert_eq!(result, child);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_child_empty() -> Result<(), Report> {
    let parent = vec![helpers::sub(b'A', 0, b'T'), helpers::sub(b'G', 5, b'C')];
    let result = compose_substitutions(&parent, &[])?;
    assert_eq!(result, parent);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_non_overlapping() -> Result<(), Report> {
    let parent = vec![helpers::sub(b'A', 0, b'T'), helpers::sub(b'A', 1, b'C')];
    let child = vec![helpers::sub(b'G', 2, b'T'), helpers::sub(b'C', 3, b'A')];
    let result = compose_substitutions(&parent, &child)?;
    let expected = vec![
      helpers::sub(b'A', 0, b'T'),
      helpers::sub(b'A', 1, b'C'),
      helpers::sub(b'G', 2, b'T'),
      helpers::sub(b'C', 3, b'A'),
    ];
    assert_eq!(result, expected);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_chain() -> Result<(), Report> {
    let parent = vec![helpers::sub(b'A', 0, b'G')];
    let child = vec![helpers::sub(b'G', 0, b'T')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(result, vec![helpers::sub(b'A', 0, b'T')]);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_cancellation() -> Result<(), Report> {
    let parent = vec![helpers::sub(b'A', 0, b'G')];
    let child = vec![helpers::sub(b'G', 0, b'A')];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(result, vec![]);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_mixed() -> Result<(), Report> {
    let parent = vec![
      helpers::sub(b'A', 0, b'T'),
      helpers::sub(b'A', 2, b'G'),
      helpers::sub(b'C', 6, b'T'),
    ];
    let child = vec![
      helpers::sub(b'G', 2, b'T'),
      helpers::sub(b'G', 4, b'C'),
      helpers::sub(b'T', 6, b'C'),
    ];
    let result = compose_substitutions(&parent, &child)?;
    let expected = vec![
      helpers::sub(b'A', 0, b'T'),
      helpers::sub(b'A', 2, b'T'),
      helpers::sub(b'G', 4, b'C'),
    ];
    assert_eq!(result, expected);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_output_sorted_by_position() -> Result<(), Report> {
    let parent = vec![
      helpers::sub(b'A', 1, b'T'),
      helpers::sub(b'G', 3, b'C'),
      helpers::sub(b'T', 5, b'A'),
    ];
    let child = vec![
      helpers::sub(b'C', 0, b'G'),
      helpers::sub(b'A', 2, b'T'),
      helpers::sub(b'C', 4, b'G'),
    ];
    let result = compose_substitutions(&parent, &child)?;
    let expected = vec![
      helpers::sub(b'C', 0, b'G'),
      helpers::sub(b'A', 1, b'T'),
      helpers::sub(b'A', 2, b'T'),
      helpers::sub(b'G', 3, b'C'),
      helpers::sub(b'C', 4, b'G'),
      helpers::sub(b'T', 5, b'A'),
    ];
    assert_eq!(result, expected);
    Ok(())
  }

  #[test]
  fn test_mutation_compose_substitutions_all_cancel() -> Result<(), Report> {
    let parent = vec![
      helpers::sub(b'A', 0, b'G'),
      helpers::sub(b'C', 1, b'T'),
      helpers::sub(b'G', 2, b'A'),
    ];
    let child = vec![
      helpers::sub(b'G', 0, b'A'),
      helpers::sub(b'T', 1, b'C'),
      helpers::sub(b'A', 2, b'G'),
    ];
    let result = compose_substitutions(&parent, &child)?;
    assert_eq!(result, vec![]);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::deletion_one_based(  b'D', PositionBase::One,  vec!["C2-", "G3-"])]
  #[case::deletion_zero_based( b'D', PositionBase::Zero, vec!["C1-", "G2-"])]
  #[case::insertion_one_based( b'I', PositionBase::One,  vec!["-2C", "-3G"])]
  #[case::insertion_zero_based(b'I', PositionBase::Zero, vec!["-1C", "-2G"])]
  #[trace]
  fn test_mutation_event_strings_expand_aligned_indels(
    #[case] kind: u8,
    #[case] positions: PositionBase,
    #[case] expected: Vec<&str>,
  ) -> Result<(), Report> {
    let segment = AlignedMutation::new((1, 3), seq![helpers::c(b'C'), helpers::c(b'G')])?;
    let event = if kind == b'D' { MutationEvent::Deletion(segment) } else { MutationEvent::Insertion(segment) };
    let actual = mutation_event_strings(&event, positions)?;
    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_based( PositionBase::One,  "A10G")]
  #[case::zero_based(PositionBase::Zero, "A9G")]
  #[trace]
  fn test_mutation_event_strings_number_substitutions(#[case] positions: PositionBase, #[case] expected: &str) -> Result<(), Report> {
    let event = MutationEvent::Substitution(helpers::sub(b'A', 9, b'G'));
    let actual = mutation_event_strings(&event, positions)?;
    assert_eq!(vec![expected], actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::canonical_change(   b'C', b'T', true)]
  #[case::incompatible_code(  b'C', b'K', true)]
  #[case::compatible_code(    b'G', b'K', false)]
  #[case::unknown(            b'C', b'N', false)]
  #[case::overlapping_codes(  b'R', b'K', false)]
  #[trace]
  fn test_mutation_sub_changes_state_only_when_state_sets_are_disjoint(
    #[case] reff: u8,
    #[case] qry: u8,
    #[case] expected: bool,
  ) -> Result<(), Report> {
    assert_eq!(expected, helpers::sub(reff, 0, qry).changes_state(&Alphabet::default())?);
    Ok(())
  }

  #[test]
  fn test_mutation_sub_changes_state_rejects_a_character_outside_the_alphabet() {
    assert_error!(
      helpers::sub(b'A', 0, b'J').changes_state(&Alphabet::default()),
      "Substitution 'A1J' holds a character outside the alphabet: Character 'J' is not in the alphabet. This is an internal error. Please report it to developers."
    );
  }

  #[test]
  fn test_mutation_edge_state_change_counts_count_only_disjoint_substitutions() -> Result<(), Report> {
    let mutation = |event| Mutation {
      track: MutationTrack::Nucleotide,
      event,
    };
    let substitution = |reff, pos, qry| mutation(MutationEvent::Substitution(helpers::sub(reff, pos, qry)));
    let deletion = mutation(MutationEvent::Deletion(AlignedMutation::new(
      (4, 5),
      seq![helpers::c(b'A')],
    )?));
    let mut graph = Graph::new();
    let root = graph.add_node();
    let [e0, e1, e2] = [(); 3].map(|()| {
      let child = graph.add_node();
      graph.add_edge(root, child).unwrap()
    });
    let edge_mutations = keymap! {
      e0 => vec![
        substitution(b'C', 0, b'T'),
        substitution(b'G', 1, b'K'),
        substitution(b'C', 2, b'N'),
        substitution(b'N', 3, b'A'),
        deletion,
      ],
      e1 => vec![substitution(b'C', 0, b'K'), substitution(b'A', 1, b'T')],
      e2 => vec![],
    };

    let actual = edge_state_change_counts(&edge_mutations, &Alphabet::default())?;

    assert_eq!(keymap! { e0 => 1, e1 => 2, e2 => 0 }, actual);
    Ok(())
  }

  proptest! {
    #[test]
    fn test_prop_mutation_event_strings_preserve_range(start in 0_usize..10_000, length in 1_usize..32) {
      let sequence = seq![helpers::c(b'A'); length];
      let event = MutationEvent::Insertion(AlignedMutation::new((start, start + length), sequence).unwrap());
      let actual = mutation_event_strings(&event, PositionBase::One).unwrap();
      let expected = (start + 1..=start + length)
        .map(|position| format!("-{position}A"))
        .collect::<Vec<_>>();
      prop_assert_eq!(expected, actual);
    }
  }

  mod helpers {
    use super::*;

    pub(super) fn c(b: u8) -> AsciiChar {
      AsciiChar::from_byte_unchecked(b)
    }

    pub(super) fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
      Sub::new(c(reff), pos, c(qry)).unwrap()
    }
  }
}
