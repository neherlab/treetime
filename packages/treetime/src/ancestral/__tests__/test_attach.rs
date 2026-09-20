use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::attach::{complete_alignment_for_leaves, sanitize_to_alphabet};
use pretty_assertions::assert_eq;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::nwk_read_str;
use treetime_primitives::{AlignmentRecord, Seq};

#[test]
fn test_attach_synthesizes_all_unknown_for_missing_tip() {
  let (graph, names) = helpers::four_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT"), ("B", "ACGT"), ("C", "ACGT")]);

  let completed = complete_alignment_for_leaves(&graph, sequences, &alphabet, false, &names).unwrap();
  let by_name = helpers::by_name(completed);

  assert_eq!(4, by_name.len());
  assert_eq!(Seq::try_from_str("ACGT").unwrap(), by_name["A"]);
  assert_eq!(Seq::try_from_str("NNNN").unwrap(), by_name["D"]);
}

#[test]
fn test_attach_aborts_above_one_third_missing() {
  let (graph, names) = helpers::two_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT")]);

  let err = complete_alignment_for_leaves(&graph, sequences, &alphabet, false, &names).unwrap_err();

  assert!(err.to_string().contains("one third"));
  assert!(err.to_string().contains("--ignore-missing-alns"));
}

#[test]
fn test_attach_ignore_missing_alns_bypasses_threshold() {
  let (graph, names) = helpers::two_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT")]);

  let completed = complete_alignment_for_leaves(&graph, sequences, &alphabet, true, &names).unwrap();
  let by_name = helpers::by_name(completed);

  assert_eq!(2, by_name.len());
  assert_eq!(Seq::try_from_str("NNNN").unwrap(), by_name["B"]);
}

#[test]
fn test_attach_exactly_one_third_missing_does_not_abort() {
  let (graph, names) = helpers::three_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT"), ("B", "ACGT")]);

  let completed = complete_alignment_for_leaves(&graph, sequences, &alphabet, false, &names).unwrap();
  let by_name = helpers::by_name(completed);

  assert_eq!(3, by_name.len());
  assert_eq!(Seq::try_from_str("NNNN").unwrap(), by_name["C"]);
}

#[test]
fn test_attach_keeps_extra_records_not_matching_any_leaf() {
  let (graph, names) = helpers::two_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT"), ("B", "ACGT"), ("reference", "ACGT")]);

  let completed = complete_alignment_for_leaves(&graph, sequences, &alphabet, false, &names).unwrap();
  let by_name = helpers::by_name(completed);

  assert_eq!(3, by_name.len());
  assert!(by_name.contains_key("reference"));
}

#[test]
fn test_sanitize_to_alphabet_folds_stop_into_unknown_for_no_stop_alphabet() {
  let aa = Alphabet::new(AlphabetName::Aa).unwrap();
  let aa_no_stop = Alphabet::new(AlphabetName::AaNoStop).unwrap();
  let seq = Seq::try_from_str("MC*X-").unwrap();

  let (kept, changed_aa) = sanitize_to_alphabet(&seq, &aa);
  assert_eq!(0, changed_aa);
  assert_eq!(seq, kept);

  let (folded, changed) = sanitize_to_alphabet(&seq, &aa_no_stop);
  assert_eq!(1, changed);
  assert_eq!(Seq::try_from_str("MCXX-").unwrap(), folded);
}

mod helpers {
  use super::*;

  pub fn two_leaf_tree() -> (Graph, BTreeMap<GraphNodeKey, Option<String>>) {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1)root;").unwrap();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    (graph, names)
  }

  pub fn three_leaf_tree() -> (Graph, BTreeMap<GraphNodeKey, Option<String>>) {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1)root;").unwrap();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    (graph, names)
  }

  pub fn four_leaf_tree() -> (Graph, BTreeMap<GraphNodeKey, Option<String>>) {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1)root;").unwrap();
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    (graph, names)
  }

  pub fn records(entries: &[(&str, &str)]) -> Vec<AlignmentRecord> {
    entries
      .iter()
      .map(|(name, seq)| AlignmentRecord {
        name: (*name).to_owned(),
        seq: Seq::try_from_str(seq).unwrap(),
      })
      .collect()
  }

  pub fn by_name(records: Vec<AlignmentRecord>) -> BTreeMap<String, Seq> {
    records.into_iter().map(|record| (record.name, record.seq)).collect()
  }
}
