#[cfg(test)]
mod tests {
  use treetime_io::nwk::nwk_fasta_node_inputs;
  use crate::commands::ancestral::aa_node_data::*;
  use crate::seq::mutation::Sub;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_utils::o;

  #[test]
  fn test_validate_aa_args_requires_cds_placeholder() {
    let err = validate_aa_args(&Some("translations.fasta".to_owned()), &["S".to_owned()], &None, &None).unwrap_err();
    assert!(err.to_string().contains("CDS placeholder"));
  }

  #[test]
  fn test_validate_aa_args_accepts_percent_gene_placeholder() {
    let result = validate_aa_args(&Some("out/%GENE.fasta".to_owned()), &["S".to_owned()], &None, &None);
    result.unwrap();
  }

  #[test]
  fn test_validate_aa_args_accepts_cds_placeholder() {
    let template = ["{", "cds", "}"].concat();
    let result = validate_aa_args(&Some(format!("out/{template}.fasta")), &["S".to_owned()], &None, &None);
    result.unwrap();
  }

  #[test]
  fn test_validate_aa_args_empty_cdses_with_annotation_ok() {
    let template = ["{", "cds", "}"].concat();
    let result = validate_aa_args(
      &Some(format!("out/{template}.fasta")),
      &[],
      &Some(PathBuf::from(concat!(env!("CARGO_MANIFEST_DIR"), "/Cargo.toml"))),
      &None,
    );
    result.unwrap();
  }

  #[test]
  fn test_validate_aa_args_empty_cdses_no_annotation_errors() {
    let template = ["{", "cds", "}"].concat();
    let err = validate_aa_args(&Some(format!("out/{template}.fasta")), &[], &None, &None).unwrap_err();
    assert!(err.to_string().contains("--cdses"));
  }

  #[allow(clippy::literal_string_with_formatting_args)]
  #[rustfmt::skip]
  #[rstest]
  #[case::cds_placeholder( "out/{cds}.fasta",  "S",  "out/S.fasta")]
  #[case::gene_placeholder("out/%GENE.fasta",   "S",  "out/S.fasta")]
  #[allow(clippy::literal_string_with_formatting_args)]
  #[case::both_placeholders("out/{cds}/%GENE.fasta", "ORF1a", "out/ORF1a/ORF1a.fasta")]
  fn test_translation_path_expands_placeholders(
    #[case] template: &str,
    #[case] cds: &str,
    #[case] expected: &str,
  ) {
    assert_eq!(PathBuf::from(expected), translation_path(template, cds));
  }

  #[test]
  fn test_validate_aa_root_sequence_cdses_requires_every_cds() {
    let by_cds = btreemap! {
      o!("S") => Seq::try_from_str("AC").unwrap(),
    };

    let err = validate_aa_root_sequence_cdses(Path::new("roots.fasta"), &by_cds, &[o!("S"), o!("M")]).unwrap_err();

    assert!(err.to_string().contains("CDS 'M'"));
  }

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
    // Both leaves carry the same sequence, so Fitch places `AC` at every node and every edge carries no
    // substitutions. The inferred root `AC` differs from the supplied reference `AA`, so only the root
    // gets a mutation (the reference-to-root diff); the leaf edges stay empty.
    let partition = helpers::fitch_partition(&graph, &names, &["AC", "AC"]);
    let reference = Seq::try_from_str("AA").unwrap();

    let actual = collect_aa_cds_node_data(&graph, &partition, "S", &names, Some(&reference)).unwrap();

    let expected = AaCdsNodeData {
      reference: o!("AA"),
      root_sequence: o!("AC"),
      node_muts: btreemap! {
        name_to_key["A"] => vec![],
        name_to_key["B"] => vec![],
        name_to_key["root"] => vec![o!("A2C")],
      },
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

  mod helpers {
    use super::*;
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::fitch::create_fitch_partition;
    use crate::ancestral::pipeline::AncestralPartition;
    use std::fmt::Write;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::fasta::read_many_fasta_str;
    use treetime_io::nwk::nwk_read_str;

    pub fn node_name_to_key(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
    ) -> BTreeMap<String, GraphNodeKey> {
      graph
        .get_nodes()
        .into_iter()
        .map(|node| {
          let node = node.read_arc();
          let key = node.key();
          let name = names[&node.key()].clone().unwrap();
          (name, key)
        })
        .collect()
    }

    pub fn named_tree() -> (Graph, BTreeMap<GraphNodeKey, Option<String>>) {
      let nwk_parsed = nwk_read_str("(A:0.1,B:0.1)root;").unwrap();
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      (graph, names)
    }

    /// Build a Fitch reconstruction whose leaves carry the given sequences (in the leaf order of the
    /// tree), for exercising the amino-acid node-data gather over a concrete partition.
    pub fn fitch_partition(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      leaf_sequences: &[&str],
    ) -> AncestralPartition {
      let alphabet = Alphabet::default();
      let mut fasta = String::new();
      for (leaf, seq) in graph.get_leaves().into_iter().zip(leaf_sequences) {
        let name = names[&leaf.read_arc().key()].clone().unwrap();
        writeln!(fasta, ">{name}").unwrap();
        writeln!(fasta, "{seq}").unwrap();
      }
      let sequences = read_many_fasta_str(&fasta, &alphabet).unwrap();
      let partition = create_fitch_partition(graph, 0, alphabet, &nwk_fasta_node_inputs(graph, names, sequences)).unwrap();
      AncestralPartition::Fitch(partition)
    }
  }
}
