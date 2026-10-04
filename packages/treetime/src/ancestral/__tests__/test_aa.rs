#[cfg(test)]
mod tests {
  use crate::alphabet::alphabet::{Alphabet, AlphabetName};
  use crate::ancestral::aa::{
    AaCdsNodeData, AaNodeData, annotation_cds_nuc_length, collect_aa_cds_node_data, diff_sequences, reconstruct_aa,
  };
  use crate::cancel::NoopCancel;
  use crate::progress::NoopProgress;
  use crate::seq::mutation::{MutationEvent, MutationTrack, Sub};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_io::nwk::nwk_read_str;
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::{assert_error, o};
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

    let mutations = partition
      .stream_sequences(&graph, &MutationTrack::AminoAcid(o!("S")), false, true, None)
      .unwrap();
    let unknown = partition.alphabet().unknown();

    let actual = collect_aa_cds_node_data(&graph, mutations, unknown, "S", Some(&reference)).unwrap();

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
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1)root;").unwrap();
    let names = nwk_parsed.names();
    let name_to_key = helpers::node_name_to_key(&names, &nwk_parsed.graph);
    let aa = Alphabet::new(AlphabetName::Aa).unwrap();
    let cdses = vec![
      helpers::cds_input("S", &aa, &[("A", "MC*"), ("B", "MC*"), ("C", "MA*")]),
      helpers::cds_input("N", &aa, &[("A", "KL"), ("B", "KL"), ("C", "KM")]),
    ];

    let actual = reconstruct_aa(
      &nwk_parsed.graph,
      &names,
      &nwk_parsed.branch_lengths,
      &helpers::sparse_params(),
      cdses,
      None,
      &NoopCancel,
      &NoopProgress,
    )
    .unwrap();

    let expected = AaNodeData {
      reference: btreemap! {
        o!("N") => o!("KL"),
        o!("S") => o!("MC*"),
      },
      root_aa_sequences: btreemap! {
        o!("N") => o!("KL"),
        o!("S") => o!("MC*"),
      },
      node_aa_mutations: btreemap! {
        name_to_key["root"] => btreemap! { o!("N") => vec![], o!("S") => vec![] },
        name_to_key["A"] => btreemap! { o!("N") => vec![], o!("S") => vec![] },
        name_to_key["B"] => btreemap! { o!("N") => vec![], o!("S") => vec![] },
        name_to_key["C"] => btreemap! {
          o!("N") => vec![helpers::substitution(b'L', 1, b'M')],
          o!("S") => vec![helpers::substitution(b'C', 1, b'A')],
        },
      },
    };
    assert_eq!(expected, actual);
    assert_eq!(vec!["N", "S"], actual.root_aa_sequences.keys().collect::<Vec<_>>());
  }

  #[test]
  fn test_reconstruct_aa_cancelled_before_first_cds_emits_nothing() {
    let (cancel, sink, emitted) = helpers::cancel_and_recording_sink(true);

    let result = helpers::reconstruct_two_cdses(&cancel, sink);

    assert_error!(result, "Operation cancelled");
    assert_eq!(Vec::<String>::new(), *emitted.lock());
  }

  #[test]
  fn test_reconstruct_aa_cancelled_after_first_cds_stops_before_second_cds() {
    let (cancel, sink, emitted) = helpers::cancel_and_recording_sink(false);

    let result = helpers::reconstruct_two_cdses(&cancel, sink);

    assert_error!(result, "Operation cancelled");
    assert_eq!(vec![o!("S"); 4], *emitted.lock());
  }

  mod helpers {
    use crate::alphabet::alphabet::{Alphabet, AlphabetName};
    use crate::ancestral::aa::{AaNodeData, AaParams, CdsInput, reconstruct_aa};
    use crate::ancestral::partition::AncestralPartition;
    use crate::cancel::Cancel;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::partition::fitch::passes::create_fitch_partition;
    use crate::partition::marginal::sample::SampleMode;
    use crate::progress::NoopProgress;
    use crate::seq::alignment::node_seq_inputs;
    use crate::seq::mutation::{MutationEvent, Sub};
    use crate::seq::sink::{SeqItem, SeqSink, SeqTrack};
    use eyre::Report;
    use parking_lot::Mutex;
    use std::collections::BTreeMap;
    use std::fmt::Write;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicBool, Ordering};
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::fasta::read_many_fasta_str;
    use treetime_io::nwk::nwk_read_str;
    use treetime_primitives::{AlignmentRecord, AsciiChar, Seq};

    pub(super) struct FlagCancel(Arc<AtomicBool>);

    impl Cancel for FlagCancel {
      fn is_cancelled(&self) -> bool {
        self.0.load(Ordering::SeqCst)
      }
    }

    struct CancellingRecordingSink {
      cancel: Arc<AtomicBool>,
      emitted: Arc<Mutex<Vec<String>>>,
    }

    impl SeqSink for CancellingRecordingSink {
      fn emit(&mut self, item: SeqItem<'_>) -> Result<(), Report> {
        let SeqTrack::Aa(name) = item.track else {
          panic!("reconstruct_aa must emit only amino acid tracks");
        };
        self.emitted.lock().push(name.to_owned());
        self.cancel.store(true, Ordering::SeqCst);
        Ok(())
      }
    }

    pub(super) fn cancel_and_recording_sink(
      initially_cancelled: bool,
    ) -> (FlagCancel, Box<dyn SeqSink>, Arc<Mutex<Vec<String>>>) {
      let flag = Arc::new(AtomicBool::new(initially_cancelled));
      let emitted = Arc::new(Mutex::new(Vec::new()));
      let sink = CancellingRecordingSink {
        cancel: Arc::clone(&flag),
        emitted: Arc::clone(&emitted),
      };
      (FlagCancel(flag), Box::new(sink), emitted)
    }

    pub(super) fn reconstruct_two_cdses(cancel: &FlagCancel, mut sink: Box<dyn SeqSink>) -> Result<AaNodeData, Report> {
      let nwk_parsed = nwk_read_str("(A:0.1,B:0.1,C:0.1)root;")?;
      let names = nwk_parsed.names();
      let aa = Alphabet::new(AlphabetName::Aa)?;
      let cdses = vec![
        cds_input("S", &aa, &[("A", "MC*"), ("B", "MC*"), ("C", "MA*")]),
        cds_input("N", &aa, &[("A", "KL"), ("B", "KL"), ("C", "KM")]),
      ];
      reconstruct_aa(
        &nwk_parsed.graph,
        &names,
        &nwk_parsed.branch_lengths,
        &sparse_params(),
        cdses,
        Some(sink.as_mut()),
        cancel,
        &NoopProgress,
      )
    }

    pub(super) fn sparse_params() -> AaParams {
      AaParams {
        dense: Some(false),
        include_leaves: false,
        impute_missing_data: false,
        sample_from_profile: SampleMode::default(),
        seed: 0,
        ignore_missing_alns: false,
      }
    }

    pub(super) fn substitution(reff: u8, pos: usize, qry: u8) -> MutationEvent {
      MutationEvent::Substitution(
        Sub::new(
          AsciiChar::from_byte_unchecked(reff),
          pos,
          AsciiChar::from_byte_unchecked(qry),
        )
        .unwrap(),
      )
    }

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
