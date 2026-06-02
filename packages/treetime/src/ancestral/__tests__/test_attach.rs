use crate::alphabet::alphabet::{Alphabet, AlphabetName};
use crate::ancestral::attach::complete_alignment_for_leaves;
use crate::payload::ancestral::GraphAncestral;
use pretty_assertions::assert_eq;
use std::collections::BTreeMap;
use treetime_io::fasta::FastaRecord;
use treetime_io::nwk::nwk_read_str;
use treetime_primitives::Seq;

#[test]
fn test_attach_synthesizes_all_unknown_for_missing_tip() {
  let graph = helpers::four_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT"), ("B", "ACGT"), ("C", "ACGT")]);

  let completed = complete_alignment_for_leaves(&graph, sequences, &alphabet, false).unwrap();
  let by_name = helpers::by_name(completed);

  assert_eq!(4, by_name.len());
  assert_eq!(Seq::try_from_str("ACGT").unwrap(), by_name["A"]);
  assert_eq!(Seq::try_from_str("NNNN").unwrap(), by_name["D"]);
}

#[test]
fn test_attach_aborts_above_one_third_missing() {
  let graph = helpers::two_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT")]);

  let err = complete_alignment_for_leaves(&graph, sequences, &alphabet, false).unwrap_err();

  assert!(err.to_string().contains("one third"));
  assert!(err.to_string().contains("--ignore-missing-alns"));
}

#[test]
fn test_attach_ignore_missing_alns_bypasses_threshold() {
  let graph = helpers::two_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT")]);

  let completed = complete_alignment_for_leaves(&graph, sequences, &alphabet, true).unwrap();
  let by_name = helpers::by_name(completed);

  assert_eq!(2, by_name.len());
  assert_eq!(Seq::try_from_str("NNNN").unwrap(), by_name["B"]);
}

#[test]
fn test_attach_exactly_one_third_missing_does_not_abort() {
  // v0 uses strict `>` with float division: 1 of 3 missing (1 > 3/3 == 1.0 is false) must not abort.
  let graph = helpers::three_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT"), ("B", "ACGT")]);

  let completed = complete_alignment_for_leaves(&graph, sequences, &alphabet, false).unwrap();
  let by_name = helpers::by_name(completed);

  assert_eq!(3, by_name.len());
  assert_eq!(Seq::try_from_str("NNNN").unwrap(), by_name["C"]);
}

#[test]
fn test_attach_keeps_extra_records_not_matching_any_leaf() {
  let graph = helpers::two_leaf_tree();
  let alphabet = Alphabet::new(AlphabetName::Nuc).unwrap();
  let sequences = helpers::records(&[("A", "ACGT"), ("B", "ACGT"), ("reference", "ACGT")]);

  let completed = complete_alignment_for_leaves(&graph, sequences, &alphabet, false).unwrap();
  let by_name = helpers::by_name(completed);

  assert_eq!(3, by_name.len());
  assert!(by_name.contains_key("reference"));
}

mod helpers {
  use super::*;

  pub fn two_leaf_tree() -> GraphAncestral {
    nwk_read_str("(A:0.1,B:0.1)root;").unwrap()
  }

  pub fn three_leaf_tree() -> GraphAncestral {
    nwk_read_str("(A:0.1,B:0.1,C:0.1)root;").unwrap()
  }

  pub fn four_leaf_tree() -> GraphAncestral {
    nwk_read_str("((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1)root;").unwrap()
  }

  pub fn records(entries: &[(&str, &str)]) -> Vec<FastaRecord> {
    entries
      .iter()
      .enumerate()
      .map(|(index, (name, seq))| FastaRecord {
        seq_name: (*name).to_owned(),
        desc: None,
        seq: Seq::try_from_str(seq).unwrap(),
        index,
      })
      .collect()
  }

  pub fn by_name(records: Vec<FastaRecord>) -> BTreeMap<String, Seq> {
    records
      .into_iter()
      .map(|record| (record.seq_name, record.seq))
      .collect()
  }
}
