#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, NON_CHAR, VARIABLE_CHAR};
  use crate::ancestral::fitch_sub::{
    finalize_sequence_forward, resolve_informative_positions_backward, resolve_nonroot_substitutions_forward,
    resolve_root_forward,
  };
  use crate::partition::sparse::{FitchSeqDistribution, SparseEdgePartition, SparseSeqInfo};
  use crate::seq::composition::Composition;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use std::sync::LazyLock;
  use treetime_primitives::{AlphabetLike, AsciiChar, BitSet128, Seq, seq, stateset};

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  fn make_seq_info(sequence: &str) -> SparseSeqInfo {
    let seq = Seq::try_from_str(sequence).unwrap();
    SparseSeqInfo {
      unknown: vec![],
      gaps: vec![],
      non_char: vec![],
      composition: Composition::new(NUC_ALPHABET.chars(), NUC_ALPHABET.gap()),
      sequence: seq,
      fitch: FitchSeqDistribution {
        variable: btreemap! {},
        variable_indel: BTreeSet::new(),
        chosen_state: btreemap! {},
      },
    }
  }

  fn make_edge() -> SparseEdgePartition {
    SparseEdgePartition::default()
  }

  #[test]
  fn test_fitch_sub_informative_backward_fixed_children_agree() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("ACGT");
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let variable = resolve_informative_positions_backward(&children, &[0], &mut sequence);

    assert!(variable.is_empty());
    assert_eq!(sequence.as_str(), "ACGT");
  }

  #[test]
  fn test_fitch_sub_informative_backward_fixed_children_disagree() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("GCGT");
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let variable = resolve_informative_positions_backward(&children, &[0], &mut sequence);

    let expected = btreemap! { 0_usize => stateset! {b'A', b'G'} };
    assert_eq!(expected, variable);
    assert_eq!(sequence[0], VARIABLE_CHAR);
  }

  #[test]
  fn test_fitch_sub_informative_backward_singleton_intersection() {
    let mut child0 = make_seq_info("ACGT");
    child0.fitch.variable.insert(0, stateset! {b'A', b'G'});
    let mut child1 = make_seq_info("ACGT");
    child1.fitch.variable.insert(0, stateset! {b'A', b'C'});
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let variable = resolve_informative_positions_backward(&children, &[0], &mut sequence);

    assert!(variable.is_empty());
    assert_eq!(sequence[0], AsciiChar::from_byte_unchecked(b'A'));
  }

  #[test]
  fn test_fitch_sub_informative_backward_multistate_intersection() {
    let mut child0 = make_seq_info("ACGT");
    child0.fitch.variable.insert(0, stateset! {b'A', b'G', b'C'});
    let mut child1 = make_seq_info("ACGT");
    child1.fitch.variable.insert(0, stateset! {b'A', b'C'});
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let variable = resolve_informative_positions_backward(&children, &[0], &mut sequence);

    let expected = btreemap! { 0_usize => stateset! {b'A', b'C'} };
    assert_eq!(expected, variable);
    assert_eq!(sequence[0], VARIABLE_CHAR);
  }

  #[test]
  fn test_fitch_sub_informative_backward_empty_intersection_uses_union() {
    let mut child0 = make_seq_info("ACGT");
    child0.fitch.variable.insert(0, stateset! {b'A'});
    let mut child1 = make_seq_info("GCGT");
    child1.fitch.variable.insert(0, stateset! {b'G'});
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let variable = resolve_informative_positions_backward(&children, &[0], &mut sequence);

    let expected = btreemap! { 0_usize => stateset! {b'A', b'G'} };
    assert_eq!(expected, variable);
    assert_eq!(sequence[0], VARIABLE_CHAR);
  }

  #[test]
  fn test_fitch_sub_informative_backward_ignores_non_char_child() {
    let mut child0 = make_seq_info("ACGT");
    child0.fitch.variable.insert(0, stateset! {b'A', b'G'});
    let mut child1 = make_seq_info("NCGT");
    child1.fitch.variable.insert(0, stateset! {b'C'});
    child1.non_char = vec![(0, 1)];
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let variable = resolve_informative_positions_backward(&children, &[0], &mut sequence);

    let expected = btreemap! { 0_usize => stateset! {b'A', b'G'} };
    assert_eq!(expected, variable);
  }

  #[test]
  fn test_fitch_sub_informative_backward_skips_local_non_char_position() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("GCGT");
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![
      NON_CHAR,
      AsciiChar::from_byte_unchecked(b'C'),
      AsciiChar::from_byte_unchecked(b'G'),
      AsciiChar::from_byte_unchecked(b'T')
    ];
    let variable = resolve_informative_positions_backward(&children, &[0], &mut sequence);

    assert!(variable.is_empty());
    assert_eq!(sequence[0], NON_CHAR);
  }

  #[test]
  fn test_fitch_sub_informative_backward_filters_transmission_child() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("GCGT");
    let mut edge0 = make_edge();
    edge0.transmission = Some(vec![(0, 1)]);
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let variable = resolve_informative_positions_backward(&children, &[0], &mut sequence);

    assert!(variable.is_empty());
    assert_eq!(sequence[0], AsciiChar::from_byte_unchecked(b'G'));
  }

  // --- resolve_root_forward ---

  #[test]
  fn test_fitch_sub_root_forward_resolves_variable() {
    let mut sequence = Seq::try_from_str("~CGT").unwrap();
    let variable = btreemap! { 0_usize => stateset! {b'A', b'G'} };
    let mut chosen_state = BTreeMap::new();

    resolve_root_forward(&mut sequence, &variable, &mut chosen_state);

    assert_eq!(
      sequence[0],
      AsciiChar::from_byte_unchecked(b'A'),
      "get_one picks alphabetically first"
    );
    assert_eq!(chosen_state[&0], AsciiChar::from_byte_unchecked(b'A'));
  }

  #[test]
  fn test_fitch_sub_root_forward_variable_indel_not_resolved() {
    // Variable indels at the root default to present (no gap).
    // Direction is resolved in the forward pass on children via parent state.
    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let variable = BTreeMap::new();
    let mut chosen_state = BTreeMap::new();

    resolve_root_forward(&mut sequence, &variable, &mut chosen_state);

    assert!(chosen_state.is_empty(), "No variable substitutions to resolve");
    assert_eq!(sequence, Seq::try_from_str("ACGT").unwrap(), "Sequence unchanged");
  }

  // --- resolve_nonroot_substitutions_forward ---

  #[test]
  fn test_fitch_sub_nonroot_forward_parent_in_child_set() -> Result<(), Report> {
    let mut sequence = Seq::try_from_str("~CGT").unwrap();
    let mut variable = btreemap! { 0_usize => stateset! {b'A', b'G'} };
    let mut chosen_state = BTreeMap::new();
    let mut composition = Composition::new(NUC_ALPHABET.chars(), NUC_ALPHABET.gap());

    let parent = make_seq_info("ACGT");

    let subs = resolve_nonroot_substitutions_forward(
      &mut sequence,
      &[],
      &mut variable,
      &mut chosen_state,
      &mut composition,
      &parent,
      &NUC_ALPHABET,
    )?;

    assert!(subs.is_empty(), "Parent state A is in child set, no substitution");
    assert_eq!(sequence[0], AsciiChar::from_byte_unchecked(b'A'));
    assert_eq!(chosen_state[&0], AsciiChar::from_byte_unchecked(b'A'));
    Ok(())
  }

  #[test]
  fn test_fitch_sub_nonroot_forward_parent_not_in_child_set() -> Result<(), Report> {
    let mut sequence = Seq::try_from_str("~CGT").unwrap();
    let mut variable = btreemap! { 0_usize => stateset! {b'C', b'G'} };
    let mut chosen_state = BTreeMap::new();
    let mut composition = Composition::new(NUC_ALPHABET.chars(), NUC_ALPHABET.gap());

    let parent = make_seq_info("ACGT");

    let subs = resolve_nonroot_substitutions_forward(
      &mut sequence,
      &[],
      &mut variable,
      &mut chosen_state,
      &mut composition,
      &parent,
      &NUC_ALPHABET,
    )?;

    assert_eq!(subs.len(), 1);
    assert_eq!(subs[0].to_string(), "A1C");
    assert_eq!(sequence[0], AsciiChar::from_byte_unchecked(b'C'));
    Ok(())
  }

  #[test]
  fn test_fitch_sub_nonroot_forward_parent_only_variable_mutation() -> Result<(), Report> {
    let mut sequence = Seq::try_from_str("GCGT").unwrap();
    let mut variable = BTreeMap::new();
    let mut chosen_state = BTreeMap::new();
    let mut composition = Composition::new(NUC_ALPHABET.chars(), NUC_ALPHABET.gap());

    let mut parent = make_seq_info("ACGT");
    parent.fitch.variable.insert(0, stateset! {b'A', b'G'});

    let subs = resolve_nonroot_substitutions_forward(
      &mut sequence,
      &[],
      &mut variable,
      &mut chosen_state,
      &mut composition,
      &parent,
      &NUC_ALPHABET,
    )?;

    assert_eq!(subs.len(), 1);
    assert_eq!(subs[0].to_string(), "A1G");
    Ok(())
  }

  #[test]
  fn test_fitch_sub_nonroot_forward_subs_sorted_by_position() -> Result<(), Report> {
    let mut sequence = Seq::try_from_str("GCGA").unwrap();
    let mut variable = btreemap! { 3_usize => stateset! {b'A'} };
    let mut chosen_state = BTreeMap::new();
    let mut composition = Composition::new(NUC_ALPHABET.chars(), NUC_ALPHABET.gap());

    let mut parent = make_seq_info("ACGT");
    parent.fitch.variable.insert(0, stateset! {b'A', b'G'});

    let subs = resolve_nonroot_substitutions_forward(
      &mut sequence,
      &[],
      &mut variable,
      &mut chosen_state,
      &mut composition,
      &parent,
      &NUC_ALPHABET,
    )?;

    assert_eq!(subs.len(), 2);
    assert_eq!(subs[0].to_string(), "A1G", "Position 0 first");
    assert_eq!(subs[1].to_string(), "T4A", "Position 3 second");
    Ok(())
  }

  // --- finalize_sequence_forward ---

  #[test]
  fn test_fitch_sub_finalize_fills_gaps_and_unknown() {
    let mut sequence = Seq::try_from_str("ACGTACGT").unwrap();
    let gaps = vec![(1, 3)];
    let unknown = vec![(5, 7)];
    let mut composition = Composition::new(NUC_ALPHABET.chars(), NUC_ALPHABET.gap());

    finalize_sequence_forward(&mut sequence, &gaps, &unknown, &mut composition, &NUC_ALPHABET, false);

    assert_eq!(sequence[1], NUC_ALPHABET.gap());
    assert_eq!(sequence[2], NUC_ALPHABET.gap());
    assert_eq!(sequence[5], NUC_ALPHABET.unknown());
    assert_eq!(sequence[6], NUC_ALPHABET.unknown());
    assert_eq!(
      sequence[0],
      AsciiChar::from_byte_unchecked(b'A'),
      "Non-gap/unknown positions unchanged"
    );
  }

  #[test]
  fn test_fitch_sub_finalize_root_recomputes_composition() {
    let mut sequence = Seq::try_from_str("AACG").unwrap();
    let mut composition = Composition::new(NUC_ALPHABET.chars(), NUC_ALPHABET.gap());

    finalize_sequence_forward(&mut sequence, &[], &[], &mut composition, &NUC_ALPHABET, true);

    let counts = composition.counts();
    assert_eq!(counts[&AsciiChar::from_byte_unchecked(b'A')], 2);
    assert_eq!(counts[&AsciiChar::from_byte_unchecked(b'C')], 1);
    assert_eq!(counts[&AsciiChar::from_byte_unchecked(b'G')], 1);
  }

  #[test]
  fn test_fitch_sub_finalize_nonroot_preserves_composition() {
    let mut sequence = Seq::try_from_str("AACG").unwrap();
    let mut composition = Composition::new(NUC_ALPHABET.chars(), NUC_ALPHABET.gap());

    finalize_sequence_forward(&mut sequence, &[], &[], &mut composition, &NUC_ALPHABET, false);

    let counts = composition.counts();
    assert_eq!(
      counts[&AsciiChar::from_byte_unchecked(b'A')],
      0,
      "Non-root does not recompute composition"
    );
  }
}
