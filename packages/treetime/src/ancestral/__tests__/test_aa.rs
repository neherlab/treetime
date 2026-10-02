#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::aa::{
    AaCdsNodeData, AaParams, annotation_cds_nuc_length, collect_aa_cds_node_data, diff_sequences, reconstruct_aa,
  };
  use crate::ancestral::sample::SampleMode;
  use crate::cancel::NoopCancel;
  use crate::progress::NoopProgress;
  use crate::seq::mutation::{MutationEvent, Sub};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::o;
  use util_augur_node_data_json::{AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotationSegment};

  #[rustfmt::skip]
  #[rstest]
  #[case::single_span(
    AugurNodeDataJsonAnnotationEntry { start: Some(100), end: Some(400), segments: None, ..Default::default() },
    Some(301)
  )]
  #[case::segments(
    AugurNodeDataJsonAnnotationEntry {
      segments: Some(vec![
        AugurNodeDataJsonAnnotationSegment { start: 1, end: 100, other: btreemap! {} },
        AugurNodeDataJsonAnnotationSegment { start: 200, end: 300, other: btreemap! {} },
      ]),
      ..Default::default()
    },
    Some(201)
  )]
  #[case::neither(
    AugurNodeDataJsonAnnotationEntry::default(),
    None
  )]
  fn test_annotation_cds_nuc_length(
    #[case] entry: AugurNodeDataJsonAnnotationEntry,
    #[case] expected: Option<i64>,
  ) {
    assert_eq!(expected, annotation_cds_nuc_length(&entry));
  }

  #[test]
  fn test_diff_sequences_skips_gap_and_unknown_states() {
    let reference = Seq::try_from_str("ACDX-").unwrap();
    let query = Seq::try_from_str("ADQXF").unwrap();

    let actual = diff_sequences(&reference, &query, AsciiChar::from_byte_unchecked(b'X')).unwrap();

    let expected = vec![
      Sub::new(
        AsciiChar::from_byte_unchecked(b'C'),
        1_usize,
        AsciiChar::from_byte_unchecked(b'D'),
      )
      .unwrap(),
      Sub::new(
        AsciiChar::from_byte_unchecked(b'D'),
        2_usize,
        AsciiChar::from_byte_unchecked(b'Q'),
      )
      .unwrap(),
    ];
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_collect_aa_cds_node_data_keeps_inferred_root_sequence() {
    let (graph, names) = helpers::named_tree();
    let name_to_key = helpers::node_name_to_key(&names, &graph);
    let partition = helpers::fitch_partition(&graph, &names, &["AC", "AC"]);
    let reference = Seq::try_from_str("AA").unwrap();

    let actual = collect_aa_cds_node_data(&graph, &partition, "S", &names, Some(&reference)).unwrap();

    let expected = AaCdsNodeData {
      reference: o!("AA"),
      root_sequence: o!("AC"),
      node_mutations: btreemap! {
        name_to_key["A"] => vec![],
        name_to_key["B"] => vec![],
        name_to_key["root"] => vec![MutationEvent::Substitution(
          Sub::new(AsciiChar::from_byte_unchecked(b'A'), 1_usize, AsciiChar::from_byte_unchecked(b'C')).unwrap()
        )],
      },
    };
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_reconstruct_aa_reconstructs_each_cds_independently_with_stop_codon() {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1)root;").unwrap();
    let names = nwk_parsed.names();
    let aa = Alphabet::new(AlphabetName::Aa).unwrap();
    let cdses = vec![
      helpers::cds_input("S", &aa, &[("A", "MC*"), ("B", "MA*")]),
      helpers::cds_input("N", &aa, &[("A", "KL"), ("B", "KM")]),
    ];
    let params = AaParams {
      dense: Some(false),
      include_leaves: false,
      impute_missing_data: false,
      sample_from_profile: SampleMode::default(),
      seed: None,
      ignore_missing_alns: false,
    };

    let actual = reconstruct_aa(
      &nwk_parsed.graph,
      &names,
      &nwk_parsed.branch_lengths,
      &params,
      cdses,
      None,
      &NoopCancel,
      &NoopProgress,
    )
    .unwrap();

    let s_root = &actual.root_aa_sequences["S"];
    assert_eq!(3, s_root.len());
    assert_eq!(Some('*'), s_root.chars().last());
    let n_root = &actual.root_aa_sequences["N"];
    assert_eq!(2, n_root.len());
    assert_eq!(Some('K'), n_root.chars().next());
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::aa::CdsInput;
    use crate::ancestral::fitch::create_fitch_partition;
    use crate::ancestral::pipeline::AncestralPartition;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::seq::alignment::node_seq_inputs;
    use std::collections::BTreeMap;
    use std::fmt::Write;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::fasta::read_many_fasta_str;
    use treetime_io::nwk::nwk_read_str;
    use treetime_primitives::{AlignmentRecord, Seq};

    pub(super) fn node_name_to_key(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
    ) -> BTreeMap<String, GraphNodeKey> {
      graph
        .get_nodes()
        .map(|node| {
          let key = node.key();
          let name = names[&node.key()].clone().unwrap();
          (name, key)
        })
        .collect()
    }

    pub(super) fn named_tree() -> (Graph, BTreeMap<GraphNodeKey, Option<String>>) {
      let nwk_parsed = nwk_read_str("(A:0.1,B:0.1)root;").unwrap();
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      (graph, names)
    }

    pub(super) fn fitch_partition(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      leaf_sequences: &[&str],
    ) -> AncestralPartition {
      let alphabet = Alphabet::default();
      let mut fasta = String::new();
      for (leaf, seq) in graph.get_leaves().zip(leaf_sequences) {
        let name = names[&leaf.key()].clone().unwrap();
        writeln!(fasta, ">{name}").unwrap();
        writeln!(fasta, "{seq}").unwrap();
      }
      let sequences: Vec<AlignmentRecord> = read_many_fasta_str(&fasta, &alphabet)
        .unwrap()
        .into_iter()
        .map(AlignmentRecord::from)
        .collect();
      let partition = create_fitch_partition(graph, 0, alphabet, &node_seq_inputs(graph, names, sequences)).unwrap();
      AncestralPartition::Fitch(partition)
    }

    pub(super) fn cds_input(name: &str, alphabet: &Alphabet, seqs: &[(&str, &str)]) -> CdsInput {
      CdsInput {
        name: name.to_owned(),
        alphabet: alphabet.clone(),
        gtr_model: GtrModelName::Infer,
        sequences: seqs
          .iter()
          .map(|(seq_name, seq)| AlignmentRecord {
            name: (*seq_name).to_owned(),
            seq: Seq::try_from_str(seq).unwrap(),
          })
          .collect(),
        annotation: None,
        reference_override: None,
      }
    }
  }
}
