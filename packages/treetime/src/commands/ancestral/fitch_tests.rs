// TODO: sparse partition tests were removed (test_fitch_internals, test_fitch_complex_gaps,
// test_fitch_deletions, test_fitch_insertions). Consider reimplementing these scenarios
// if sparse mode coverage is insufficient.
use crate::alphabet::alphabet::Alphabet;
use crate::commands::ancestral::fitch::{ancestral_reconstruction_fitch, compress_sequences, get_common_length};
use crate::io::fasta::read_many_fasta_str;
use crate::io::nwk::nwk_read_str;
use crate::representation::graph_ancestral::GraphAncestral;
use crate::representation::partition_fitch::PartitionFitch;
use eyre::Report;
use indoc::indoc;
use lazy_static::lazy_static;
use maplit::btreemap;
use parking_lot::RwLock;
use std::collections::BTreeMap;
use std::sync::Arc;
use treetime_io::json::{JsonPretty, json_write_str};

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

  let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

  let alphabet = Alphabet::default();

  let partitions_parsimony = [Arc::new(RwLock::new(PartitionFitch {
    index: 0,
    alphabet,
    length: get_common_length(&aln)?,
    nodes: btreemap! {},
    edges: btreemap! {},
  }))];

  compress_sequences(&graph, &partitions_parsimony, &aln)?;

  let mut actual = BTreeMap::new();
  ancestral_reconstruction_fitch(&graph, false, &partitions_parsimony, |node, seq| {
    actual.insert(node.payload.name.clone(), seq.to_string());
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

  let graph: GraphAncestral = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

  let alphabet = Alphabet::default();

  let partitions_parsimony = [Arc::new(RwLock::new(PartitionFitch {
    index: 0,
    alphabet,
    length: get_common_length(&aln)?,
    nodes: btreemap! {},
    edges: btreemap! {},
  }))];

  compress_sequences(&graph, &partitions_parsimony, &aln)?;

  let mut actual = BTreeMap::new();
  ancestral_reconstruction_fitch(&graph, true, &partitions_parsimony, |node, seq| {
    actual.insert(node.payload.name.clone(), seq.to_string());
  })?;

  assert_eq!(
    json_write_str(&expected, JsonPretty(false))?,
    json_write_str(&actual, JsonPretty(false))?
  );

  Ok(())
}
