#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
  use crate::ancestral::fitch_sub::{
    finalize_sequence_forward, resolve_fixed_positions_backward, resolve_nonroot_substitutions_forward,
    resolve_root_forward, resolve_variable_positions_backward,
  };
  use crate::representation::payload::sparse::{Deletion, FitchSeqDistribution, SparseEdgePartition, SparseSeqInfo};
  use crate::seq::composition::Composition;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
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
        variable_indel: btreemap! {},
        composition: Composition::new(NUC_ALPHABET.chars(), NUC_ALPHABET.gap()),
        chosen_state: btreemap! {},
      },
    }
  }

  fn make_edge() -> SparseEdgePartition {
    SparseEdgePartition::default()
  }

  // --- resolve_variable_positions_backward ---

  #[test]
  fn test_fitch_sub_variable_backward_children_agree() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("ACGT");
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &mut sequence);

    assert!(variable.is_empty(), "No variable positions when children agree");
  }

  #[test]
  fn test_fitch_sub_variable_backward_children_disagree() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("GCGT");
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &mut sequence);

    assert!(variable.is_empty(), "Fixed-position disagreement is not handled here");
    assert_eq!(sequence[1], FILL_CHAR, "Position 1 unchanged by variable pass");
  }

  #[test]
  fn test_fitch_sub_variable_backward_intersection_ambiguous() {
    let mut child0 = make_seq_info("ACGT");
    child0.fitch.variable.insert(0, stateset! {b'A', b'G'});
    let mut child1 = make_seq_info("ACGT");
    child1.fitch.variable.insert(0, stateset! {b'A', b'C'});
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &mut sequence);

    // intersection {A,G} ∩ {A,C} = {A} is a singleton: resolved immediately, not stored as variable
    assert!(variable.is_empty());
    assert_eq!(sequence[0], AsciiChar::from_byte_unchecked(b'A'));
  }

  #[test]
  fn test_fitch_sub_variable_backward_intersection_multi_element() {
    let mut child0 = make_seq_info("ACGT");
    child0.fitch.variable.insert(0, stateset! {b'A', b'G', b'C'});
    let mut child1 = make_seq_info("ACGT");
    child1.fitch.variable.insert(0, stateset! {b'A', b'C'});
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &mut sequence);

    // intersection {A,G,C} ∩ {A,C} = {A,C}: ambiguous, stored as variable
    assert_eq!(variable.len(), 1);
    assert!(variable[&0].contains(b'A'));
    assert!(variable[&0].contains(b'C'));
    assert!(!variable[&0].contains(b'G'));
    assert_eq!(sequence[0], VARIABLE_CHAR);
  }

  #[test]
  fn test_fitch_sub_variable_backward_intersection_empty_takes_union() {
    let mut child0 = make_seq_info("ACGT");
    child0.fitch.variable.insert(0, stateset! {b'A'});
    let mut child1 = make_seq_info("GCGT");
    child1.fitch.variable.insert(0, stateset! {b'G'});
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &mut sequence);

    assert_eq!(variable.len(), 1);
    assert!(variable[&0].contains(b'A'));
    assert!(variable[&0].contains(b'G'));
    assert_eq!(sequence[0], VARIABLE_CHAR);
  }

  #[test]
  fn test_fitch_sub_variable_backward_singleton_intersection_resolved() {
    let mut child0 = make_seq_info("ACGT");
    child0.fitch.variable.insert(0, stateset! {b'A', b'G'});
    let mut child1 = make_seq_info("ACGT");
    child1.fitch.variable.insert(0, stateset! {b'A'});
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &mut sequence);

    assert!(variable.is_empty(), "Singleton intersection resolved immediately");
    assert_eq!(sequence[0], AsciiChar::from_byte_unchecked(b'A'));
  }

  #[test]
  fn test_fitch_sub_variable_backward_skips_non_char_child() {
    let mut child0 = make_seq_info("ACGT");
    child0.fitch.variable.insert(0, stateset! {b'A', b'G'});
    let mut child1 = make_seq_info("NCGT");
    child1.fitch.variable.insert(0, stateset! {b'C'});
    child1.non_char = vec![(0, 1)];
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &mut sequence);

    // child1 is in non_char at pos 0, so only child0's {A,G} contributes
    assert_eq!(variable.len(), 1);
    assert!(variable[&0].contains(b'A'));
    assert!(variable[&0].contains(b'G'));
  }

  // --- resolve_fixed_positions_backward ---

  #[test]
  fn test_fitch_sub_fixed_backward_sets_fill_char() {
    let child0 = make_seq_info("ACGT");
    let edge0 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0)];

    let mut sequence = seq![FILL_CHAR; 4];
    let mut variable = BTreeMap::new();
    resolve_fixed_positions_backward(&children, &NUC_ALPHABET, &mut sequence, &mut variable);

    assert_eq!(sequence.as_str(), "ACGT");
    assert!(variable.is_empty());
  }

  #[test]
  fn test_fitch_sub_fixed_backward_promotes_to_variable() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("GCGT");
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let mut variable = BTreeMap::new();
    resolve_fixed_positions_backward(&children, &NUC_ALPHABET, &mut sequence, &mut variable);

    assert_eq!(variable.len(), 1);
    assert!(variable[&0].contains(b'A'));
    assert!(variable[&0].contains(b'G'));
    assert_eq!(sequence[0], VARIABLE_CHAR);
    assert_eq!(sequence[1], AsciiChar::from_byte_unchecked(b'C'));
  }

  #[test]
  fn test_fitch_sub_fixed_backward_skips_non_char() {
    let child0 = make_seq_info("ACGT");
    let edge0 = make_edge();
    let children: Vec<(&SparseSeqInfo, &SparseEdgePartition)> = vec![(&child0, &edge0)];

    let mut sequence = seq![NON_CHAR; 4];
    let mut variable = BTreeMap::new();
    resolve_fixed_positions_backward(&children, &NUC_ALPHABET, &mut sequence, &mut variable);

    assert_eq!(sequence[0], NON_CHAR, "NON_CHAR positions are skipped");
    assert!(variable.is_empty());
  }

  // --- resolve_root_forward ---

  #[test]
  fn test_fitch_sub_root_forward_resolves_variable() {
    let mut sequence = Seq::try_from_str("~CGT").unwrap();
    let mut gaps = vec![];
    let variable = btreemap! { 0_usize => stateset! {b'A', b'G'} };
    let variable_indel = BTreeMap::new();
    let mut chosen_state = BTreeMap::new();

    resolve_root_forward(&mut sequence, &mut gaps, &variable, &variable_indel, &mut chosen_state);

    assert_eq!(
      sequence[0],
      AsciiChar::from_byte_unchecked(b'A'),
      "get_one picks alphabetically first"
    );
    assert_eq!(chosen_state[&0], AsciiChar::from_byte_unchecked(b'A'));
  }

  #[test]
  fn test_fitch_sub_root_forward_majority_gap() {
    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let mut gaps = vec![];
    let variable = BTreeMap::new();
    let variable_indel = btreemap! { (1_usize, 3_usize) => Deletion { deleted: 3, present: 1 } };
    let mut chosen_state = BTreeMap::new();

    resolve_root_forward(&mut sequence, &mut gaps, &variable, &variable_indel, &mut chosen_state);

    assert_eq!(gaps, vec![(1, 3)], "Majority deleted -> gap added");
  }

  #[test]
  fn test_fitch_sub_root_forward_minority_gap_not_added() {
    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let mut gaps = vec![];
    let variable = BTreeMap::new();
    let variable_indel = btreemap! { (1_usize, 3_usize) => Deletion { deleted: 1, present: 3 } };
    let mut chosen_state = BTreeMap::new();

    resolve_root_forward(&mut sequence, &mut gaps, &variable, &variable_indel, &mut chosen_state);

    assert!(gaps.is_empty(), "Minority deleted -> no gap");
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
