#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName, FILL_CHAR, NON_CHAR, VARIABLE_CHAR};
  use crate::ancestral::fitch_sub::{
    discover_fixed_disagreements_backward, finalize_sequence_forward, resolve_nonroot_substitutions_forward,
    resolve_root_forward, resolve_variable_positions_backward,
  };
  use crate::partition::storage::sparse::{FitchSeqDistribution, FitchSeqInfo, SparseEdgeObs};
  use crate::seq::composition::Composition;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::{BTreeMap, BTreeSet};
  use std::sync::LazyLock;
  use treetime_primitives::{AlphabetLike, AsciiChar, BitSet128, Seq, seq, stateset};

  static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  fn make_seq_info(sequence: &str) -> FitchSeqInfo {
    let seq = Seq::try_from_str(sequence).unwrap();
    FitchSeqInfo {
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

  fn make_edge() -> SparseEdgeObs {
    SparseEdgeObs::default()
  }

  #[test]
  fn test_fitch_sub_variable_backward_children_agree() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("ACGT");
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &[], &mut sequence);

    assert!(variable.is_empty(), "No variable positions when children agree");
  }

  #[test]
  fn test_fitch_sub_variable_backward_children_disagree() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("GCGT");
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &[], &mut sequence);

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
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &[], &mut sequence);

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
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &[], &mut sequence);

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
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &[], &mut sequence);

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
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &[], &mut sequence);

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
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let variable = resolve_variable_positions_backward(&children, &[], &[], &mut sequence);

    assert_eq!(variable.len(), 1);
    assert!(variable[&0].contains(b'A'));
    assert!(variable[&0].contains(b'G'));
  }

  #[test]
  fn test_fitch_sub_discover_sets_fill_char() {
    let child0 = make_seq_info("ACGT");
    let edge0 = make_edge();
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0)];

    let mut sequence = seq![FILL_CHAR; 4];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);

    assert_eq!(sequence.as_str(), "ACGT");
    assert!(discovered.is_empty());
  }

  #[test]
  fn test_fitch_sub_discover_flags_disagreement() {
    let child0 = make_seq_info("ACGT");
    let child1 = make_seq_info("GCGT");
    let edge0 = make_edge();
    let edge1 = make_edge();
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0), (&child1, &edge1)];

    let mut sequence = seq![FILL_CHAR; 4];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);

    assert_eq!(discovered, vec![0], "only the disagreeing position is reported");
    assert_eq!(sequence[0], VARIABLE_CHAR);
    assert_eq!(sequence[1], AsciiChar::from_byte_unchecked(b'C'));
  }

  #[test]
  fn test_fitch_sub_discover_skips_non_char() {
    let child0 = make_seq_info("ACGT");
    let edge0 = make_edge();
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![(&child0, &edge0)];

    let mut sequence = seq![NON_CHAR; 4];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);

    assert_eq!(sequence[0], NON_CHAR, "NON_CHAR positions are skipped");
    assert!(discovered.is_empty());
  }

  #[test]
  fn test_fitch_sub_discover_reports_position_once_for_many_children() {
    let child0 = make_seq_info("A");
    let child1 = make_seq_info("C");
    let child2 = make_seq_info("G");
    let child3 = make_seq_info("T");
    let edges = [make_edge(), make_edge(), make_edge(), make_edge()];
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> = vec![
      (&child0, &edges[0]),
      (&child1, &edges[1]),
      (&child2, &edges[2]),
      (&child3, &edges[3]),
    ];

    let mut sequence = seq![FILL_CHAR; 1];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);

    assert_eq!(discovered, vec![0]);
  }

  fn fixed_children(states: &[&str]) -> Vec<FitchSeqInfo> {
    states.iter().map(|s| make_seq_info(s)).collect()
  }

  fn as_children<'a>(
    infos: &'a [FitchSeqInfo],
    edges: &'a [SparseEdgeObs],
  ) -> Vec<(&'a FitchSeqInfo, &'a SparseEdgeObs)> {
    infos.iter().zip(edges.iter()).collect()
  }

  #[test]
  fn test_fitch_sub_plurality_three_children_c_c_a() {
    let infos = fixed_children(&["C", "C", "A"]);
    let edges = vec![make_edge(), make_edge(), make_edge()];
    let children = as_children(&infos, &edges);

    let mut sequence = seq![FILL_CHAR; 1];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);
    let variable = resolve_variable_positions_backward(&children, &discovered, &[], &mut sequence);

    assert_eq!(variable[&0], stateset! {b'C'}, "only C is of minimum cost");
    assert_eq!(sequence[0], VARIABLE_CHAR);
  }

  #[test]
  fn test_fitch_sub_plurality_five_children_four_to_one() {
    let infos = fixed_children(&["C", "C", "C", "C", "A"]);
    let edges = vec![make_edge(), make_edge(), make_edge(), make_edge(), make_edge()];
    let children = as_children(&infos, &edges);

    let mut sequence = seq![FILL_CHAR; 1];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);
    let variable = resolve_variable_positions_backward(&children, &discovered, &[], &mut sequence);

    assert_eq!(variable[&0], stateset! {b'C'});
  }

  #[test]
  fn test_fitch_sub_plurality_retains_genuine_tie() {
    let infos = fixed_children(&["C", "C", "A", "A"]);
    let edges = vec![make_edge(), make_edge(), make_edge(), make_edge()];
    let children = as_children(&infos, &edges);

    let mut sequence = seq![FILL_CHAR; 1];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);
    let variable = resolve_variable_positions_backward(&children, &discovered, &[], &mut sequence);

    assert_eq!(variable[&0], stateset! {b'A', b'C'}, "an even split stays ambiguous");
  }

  #[test]
  fn test_fitch_sub_plurality_two_children_unchanged() {
    let infos = fixed_children(&["C", "A"]);
    let edges = vec![make_edge(), make_edge()];
    let children = as_children(&infos, &edges);

    let mut sequence = seq![FILL_CHAR; 1];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);
    let variable = resolve_variable_positions_backward(&children, &discovered, &[], &mut sequence);

    assert_eq!(variable[&0], stateset! {b'A', b'C'});
  }

  #[test]
  fn test_fitch_sub_plurality_counts_ambiguous_leaf_exactly_once() {
    let mut child0 = make_seq_info("C");
    let mut child1 = make_seq_info("A");
    let mut child2 = make_seq_info("R");
    child2
      .fitch
      .variable
      .insert(0, NUC_ALPHABET.char_to_set(AsciiChar::from_byte_unchecked(b'R')));
    child0.fitch.variable.clear();
    child1.fitch.variable.clear();
    let edges = [make_edge(), make_edge(), make_edge()];
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> =
      vec![(&child0, &edges[0]), (&child1, &edges[1]), (&child2, &edges[2])];

    let mut sequence = seq![FILL_CHAR; 1];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);
    assert_eq!(discovered, vec![0], "the fixed disagreement is discovered");

    let variable = resolve_variable_positions_backward(&children, &discovered, &[], &mut sequence);

    assert_eq!(variable[&0], stateset! {b'A'});
  }

  #[test]
  fn test_fitch_sub_plurality_position_from_both_sources_resolved_once() {
    let mut child0 = make_seq_info("C");
    child0.fitch.variable.insert(0, stateset! {b'C'});
    let child1 = make_seq_info("A");
    let child2 = make_seq_info("A");
    let edges = [make_edge(), make_edge(), make_edge()];
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> =
      vec![(&child0, &edges[0]), (&child1, &edges[1]), (&child2, &edges[2])];

    let mut sequence = seq![FILL_CHAR; 1];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);
    let variable = resolve_variable_positions_backward(&children, &discovered, &[], &mut sequence);

    assert_eq!(variable.len(), 1, "the position appears once in the merged list");
    assert_eq!(variable[&0], stateset! {b'A'}, "count(A) = 2 beats count(C) = 1");
  }

  #[test]
  fn test_fitch_sub_resolution_never_stores_sentinel_states() {
    let mut child0 = make_seq_info("CA");
    child0.fitch.variable.insert(1, stateset! {b'A', b'G'});
    let child1 = make_seq_info("AA");
    let child2 = make_seq_info("GA");
    let edges = [make_edge(), make_edge(), make_edge()];
    let children: Vec<(&FitchSeqInfo, &SparseEdgeObs)> =
      vec![(&child0, &edges[0]), (&child1, &edges[1]), (&child2, &edges[2])];

    let mut sequence = seq![FILL_CHAR; 2];
    let discovered = discover_fixed_disagreements_backward(&children, &NUC_ALPHABET, &mut sequence);
    let variable = resolve_variable_positions_backward(&children, &discovered, &[], &mut sequence);

    for (pos, states) in &variable {
      for state in states.chars() {
        assert!(
          NUC_ALPHABET.is_canonical(state) || NUC_ALPHABET.is_ambiguous(state),
          "position {pos} holds non-alphabet state {state:?}"
        );
        assert_ne!(VARIABLE_CHAR, state);
        assert_ne!(FILL_CHAR, state);
        assert_ne!(NON_CHAR, state);
      }
    }
  }

  #[test]
  fn test_fitch_sub_root_forward_resolves_variable() {
    let mut sequence = Seq::try_from_str("~CGT").unwrap();
    let variable = btreemap! { 0_usize => stateset! {b'A', b'G'} };
    let mut chosen_state = BTreeMap::new();

    resolve_root_forward(&mut sequence, &variable, &mut chosen_state, &NUC_ALPHABET);

    assert_eq!(
      sequence[0],
      AsciiChar::from_byte_unchecked(b'A'),
      "ties break to the first canonical state"
    );
    assert_eq!(chosen_state[&0], AsciiChar::from_byte_unchecked(b'A'));
  }

  #[test]
  fn test_fitch_sub_root_forward_tie_break_uses_canonical_order_not_ascii() {
    let aa = Alphabet::new(AlphabetName::Aa).unwrap();
    let mut sequence = Seq::try_from_str("~").unwrap();
    let variable = btreemap! { 0_usize => stateset! {b'*', b'A', b'C'} };
    let mut chosen_state = BTreeMap::new();

    resolve_root_forward(&mut sequence, &variable, &mut chosen_state, &aa);

    assert_eq!(sequence[0], AsciiChar::from_byte_unchecked(b'A'));
    assert_eq!(chosen_state[&0], AsciiChar::from_byte_unchecked(b'A'));
    assert_eq!(
      variable[&0].get_one(),
      AsciiChar::from_byte_unchecked(b'*'),
      "get_one would have selected the stop codon"
    );
  }

  #[test]
  fn test_fitch_sub_root_forward_tie_break_is_deterministic() {
    let variable = btreemap! { 0_usize => stateset! {b'A', b'C', b'G', b'T'} };
    for _ in 0..16 {
      let mut sequence = Seq::try_from_str("~").unwrap();
      let mut chosen_state = BTreeMap::new();
      resolve_root_forward(&mut sequence, &variable, &mut chosen_state, &NUC_ALPHABET);
      assert_eq!(sequence[0], AsciiChar::from_byte_unchecked(b'A'));
    }
  }

  #[test]
  fn test_fitch_sub_root_forward_variable_indel_not_resolved() {
    let mut sequence = Seq::try_from_str("ACGT").unwrap();
    let variable = BTreeMap::new();
    let mut chosen_state = BTreeMap::new();

    resolve_root_forward(&mut sequence, &variable, &mut chosen_state, &NUC_ALPHABET);

    assert!(chosen_state.is_empty(), "No variable substitutions to resolve");
    assert_eq!(sequence, Seq::try_from_str("ACGT").unwrap(), "Sequence unchanged");
  }

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
