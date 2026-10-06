#[cfg(test)]
mod tests {
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_utils::io::json::json_read_str;
  use treetime_utils::o;
  use util_augur_node_data_json::AugurNodeDataJsonAncestral;

  #[test]
  fn test_augur_node_data_ancestral_full_output() {
    let (graph, names, maps) = helpers::mutation_case();
    let actual = helpers::write_json(&graph, &names, maps, &[false, false, false, false]);

    let expected = format!(
      r#"{{
  "generated_by": {{
    "program": "treetime",
    "version": "{version}"
  }},
  "nodes": {{
    "A": {{
      "muts": [
        "T4A"
      ],
      "sequence": "ACGA"
    }},
    "B": {{
      "muts": [],
      "sequence": "ACGT"
    }},
    "root": {{
      "muts": [],
      "sequence": "ACGT"
    }}
  }},
  "annotations": {{
    "nuc": {{
      "start": 1,
      "end": 4,
      "strand": "+",
      "type": "source"
    }}
  }},
  "reference": {{
    "nuc": "ACGT"
  }},
  "mask": "0000"
}}"#,
      version = env!("CARGO_PKG_VERSION")
    );

    assert_eq!(expected, actual.trim());
  }

  #[test]
  fn test_augur_node_data_ancestral_roundtrip() {
    let (graph, names, maps) = helpers::mutation_case();
    let json_str = helpers::write_json(&graph, &names, maps, &[false, false, false, false]);

    let original: serde_json::Value = serde_json::from_str(&json_str).unwrap();
    let typed: AugurNodeDataJsonAncestral = json_read_str(&json_str).unwrap();
    let roundtripped: serde_json::Value = serde_json::to_value(&typed).unwrap();

    assert_eq!(original, roundtripped);
  }

  #[test]
  fn test_augur_node_data_ancestral_mask_filters_mutations() {
    let (graph, names, maps) = helpers::mutation_case();
    let actual = helpers::write_json(&graph, &names, maps, &[false, false, false, true]);

    let expected = format!(
      r#"{{
  "generated_by": {{
    "program": "treetime",
    "version": "{version}"
  }},
  "nodes": {{
    "A": {{
      "muts": [],
      "sequence": "ACGN"
    }},
    "B": {{
      "muts": [],
      "sequence": "ACGN"
    }},
    "root": {{
      "muts": [],
      "sequence": "ACGN"
    }}
  }},
  "annotations": {{
    "nuc": {{
      "start": 1,
      "end": 4,
      "strand": "+",
      "type": "source"
    }}
  }},
  "reference": {{
    "nuc": "ACGT"
  }},
  "mask": "0001"
}}"#,
      version = env!("CARGO_PKG_VERSION")
    );

    assert_eq!(expected, actual.trim());
  }

  #[test]
  fn test_augur_node_data_ancestral_root_has_empty_muts() {
    let (graph, names, maps) = helpers::mutation_case();
    let json_str = helpers::write_json(&graph, &names, maps, &[false, false, false, false]);
    let data: AugurNodeDataJsonAncestral = json_read_str(&json_str).unwrap();

    assert_eq!(Vec::<String>::new(), data.nodes["root"].muts);
    assert_eq!(vec!["T4A".to_owned()], data.nodes["A"].muts);
  }

  #[test]
  fn test_augur_node_data_ancestral_parsimony_end_to_end() {
    use crate::commands::shared::method_anc::MethodAncestralCli;
    use crate::commands::shared::model::GtrModelNameCli;
    let actual = helpers::reconstruct_json(MethodAncestralCli::Parsimony, None, GtrModelNameCli::Infer);
    assert_eq!(helpers::expected_invariant_json(), actual.trim());
  }

  #[test]
  fn test_augur_node_data_ancestral_marginal_sparse_end_to_end() {
    use crate::commands::shared::method_anc::MethodAncestralCli;
    use crate::commands::shared::model::GtrModelNameCli;
    let actual = helpers::reconstruct_json(MethodAncestralCli::Marginal, Some(false), GtrModelNameCli::JC69);
    assert_eq!(helpers::expected_invariant_json(), actual.trim());
  }

  #[test]
  fn test_augur_node_data_ancestral_marginal_dense_end_to_end() {
    use crate::commands::shared::method_anc::MethodAncestralCli;
    use crate::commands::shared::model::GtrModelNameCli;
    let actual = helpers::reconstruct_json(MethodAncestralCli::Marginal, Some(true), GtrModelNameCli::JC69);
    assert_eq!(helpers::expected_invariant_json(), actual.trim());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::sparse(false)]
  #[case::dense( true)]
  #[trace]
  fn test_augur_node_data_ancestral_mask_ignores_records_not_in_tree(#[case] dense: bool) {
    let data = helpers::reconstruct_json_with_outlier_record(dense);

    let sequences = data
      .nodes
      .iter()
      .filter_map(|(name, node)| Some((name.clone(), node.sequence.clone()?)))
      .collect::<BTreeMap<_, _>>();
    assert_eq!(
      (
        Some(o!("001100")),
        btreemap! {
          o!("A") => o!("ACNN-T"),
          o!("B") => o!("ACNN-T"),
          o!("C") => o!("ACNN-T"),
          o!("AB") => o!("ACNN-T"),
          o!("root") => o!("ACNN-T"),
        }
      ),
      (data.metadata.mask, sequences)
    );
  }

  #[test]
  fn test_augur_node_data_ancestral_with_aa_reconstruction() {
    let actual = helpers::build_json_with_aa();
    let expected = helpers::expected_json_with_aa();

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_augur_node_data_ancestral_multi_cds_translations_end_to_end() {
    let data = helpers::reconstruct_json_with_translations();

    assert_eq!(
      Some(btreemap! { o!("M") => o!("WY"), o!("S") => o!("MKL") }),
      data.nodes["root"].aa_sequences
    );

    let no_muts = btreemap! { o!("M") => Vec::<String>::new(), o!("S") => Vec::<String>::new() };
    assert_eq!(Some(no_muts.clone()), data.nodes["A"].aa_muts);
    assert_eq!(Some(no_muts), data.nodes["B"].aa_muts);
  }

  mod helpers {
    use crate::__tests__::test_support::tests::NUC_ALPHABET;
    use crate::commands::ancestral::args::{TreetimeAncestralArgs, TreetimeAncestralArgsRaw};
    use crate::commands::ancestral::run::run_ancestral_reconstruction;
    use crate::commands::shared::alignment::AlignmentArgs;
    use crate::commands::shared::method_anc::MethodAncestralCli;
    use crate::commands::shared::model::GtrModelNameCli;
    use crate::commands::shared::model::ModelArgs;
    use crate::commands::shared::output_args::OutputCoreArgs;
    use app_output::annotated_graph::{AnnotatedGraph, AnnotatedTreeView, Divergence, TreeAminoAcids, TreeSequences};
    use app_output::augur_node_data_ancestral::{AncestralNodeSequences, build_augur_node_data_ancestral};
    use maplit::btreemap;
    use std::collections::BTreeMap;
    use tempfile::tempdir;
    use treetime::alphabet::alphabet::Alphabet;
    use treetime::ancestral::aa::{AaCdsNodeData, AaNodeData};
    use treetime::cancel::NoopCancel;
    use treetime::progress::NoopProgress;
    use treetime::seq::mutation::{Mutation, MutationEvent, MutationTrack, Sub};
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::{AsciiChar, Seq};
    use treetime_utils::io::json::{JsonPretty, json_read_str, json_write_str};
    use treetime_utils::o;
    use util_augur_node_data_json::{
      AugurNodeDataJsonAncestral, AugurNodeDataJsonAncestralMeta, AugurNodeDataJsonAncestralNode,
      AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotations, AugurNodeDataJsonGeneratedBy,
    };

    fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
      Sub::new(
        AsciiChar::from_byte_unchecked(reff),
        pos,
        AsciiChar::from_byte_unchecked(qry),
      )
      .unwrap()
    }

    pub(super) fn mutation_case() -> (Graph, BTreeMap<GraphNodeKey, Option<String>>, OutputMaps) {
      let nwk_parsed = nwk_read(b"(A:0.1,B:0.1)root;".as_slice()).unwrap();
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let graph: Graph = graph;
      let seqs = btreemap! { o!("A") => o!("ACGA"), o!("B") => o!("ACGT"), o!("root") => o!("ACGT") };
      let edge_subs = btreemap! { o!("A") => vec![sub(b'T', 3, b'A')] };
      let maps = build_output_maps(&graph, &names, &seqs, &edge_subs, 4);
      (graph, names, maps)
    }

    fn node_name_to_key(
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

    pub(super) struct OutputMaps {
      root_sequence: Seq,
      edge_mutations: BTreeMap<GraphEdgeKey, Vec<Mutation>>,
      node_sequences: BTreeMap<GraphNodeKey, Seq>,
      length: usize,
    }

    fn build_output_maps(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      seqs: &BTreeMap<String, String>,
      edge_subs_by_child: &BTreeMap<String, Vec<Sub>>,
      length: usize,
    ) -> OutputMaps {
      let name_of = |key: GraphNodeKey| names[&key].clone().unwrap();
      let node_sequences = graph
        .get_nodes()
        .map(|node| (node.key(), Seq::try_from_str(&seqs[&name_of(node.key())]).unwrap()))
        .collect();
      let edge_mutations = graph
        .get_edges()
        .map(|edge| {
          let subs = edge_subs_by_child
            .get(&name_of(edge.target()))
            .cloned()
            .unwrap_or_default();
          let mutations = subs
            .into_iter()
            .map(|sub| Mutation::substitution(MutationTrack::Nucleotide, sub))
            .collect();
          (edge.key(), mutations)
        })
        .collect();
      OutputMaps {
        root_sequence: Seq::try_from_str(&seqs["root"]).unwrap(),
        edge_mutations,
        node_sequences,
        length,
      }
    }

    fn build_json(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      maps: OutputMaps,
      mask: &[bool],
      amino_acids: Option<TreeAminoAcids<'_>>,
    ) -> AugurNodeDataJsonAncestral {
      let branch_lengths = graph.get_edges().map(|edge| (edge.key(), None)).collect();
      let annotated = AnnotatedGraph {
        graph,
        names,
        divergence_branch_lengths: &branch_lengths,
        time_branch_lengths: None,
        divergence: Divergence::CumulativeBranchLength,
        sequences: Some(TreeSequences {
          alphabet: &NUC_ALPHABET,
          root_sequence: &maps.root_sequence,
          edge_mutations: &maps.edge_mutations,
          mutation_counts: None,
          amino_acids,
        }),
        dates: None,
        traits: None,
      };
      let sequences = AncestralNodeSequences {
        node_sequences: maps.node_sequences,
        alignment_length: maps.length,
        ambiguous_char: Alphabet::default().unknown(),
        mask,
      };
      build_augur_node_data_ancestral(&AnnotatedTreeView::new(&annotated).unwrap(), sequences).unwrap()
    }

    pub(super) fn write_json(
      graph: &Graph,
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      maps: OutputMaps,
      mask: &[bool],
    ) -> String {
      let data = build_json(graph, names, maps, mask, None);
      json_write_str(&data, JsonPretty(true)).unwrap()
    }

    pub(super) fn reconstruct_json(method: MethodAncestralCli, dense: Option<bool>, model: GtrModelNameCli) -> String {
      let dir = tempdir().unwrap();
      let tree_path = dir.path().join("tree.nwk");
      let fasta_path = dir.path().join("aln.fasta");
      let node_data_path = dir.path().join("augur-node-data.json");
      std::fs::write(&tree_path, "(A:0.1,B:0.1)root;").unwrap();
      std::fs::write(&fasta_path, ">A\nACGT\n>B\nACGT\n").unwrap();

      let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![fasta_path],
        },
        tree: Some(tree_path),
        method_anc: method,
        dense,
        model_args: ModelArgs {
          model,
          ..ModelArgs::default()
        },
        output: OutputCoreArgs {
          output_tree_nwk: Some(dir.path().join("tree_out.nwk")),
          ..Default::default()
        },
        output_augur_node_data: Some(node_data_path.clone()),
        ..TreetimeAncestralArgsRaw::default()
      })
      .unwrap();

      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress).unwrap();
      std::fs::read_to_string(node_data_path).unwrap()
    }

    pub(super) fn reconstruct_json_with_outlier_record(dense: bool) -> AugurNodeDataJsonAncestral {
      let dir = tempdir().unwrap();
      let tree_path = dir.path().join("tree.nwk");
      let fasta_path = dir.path().join("aln.fasta");
      let node_data_path = dir.path().join("augur-node-data.json");
      std::fs::write(&tree_path, "((A:0.1,B:0.1)AB:0.1,C:0.1)root;").unwrap();
      std::fs::write(&fasta_path, ">A\nACNN-T\n>B\nACNN-T\n>C\nACNN-T\n>outlier\nAAGTGT\n").unwrap();

      let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![fasta_path],
        },
        tree: Some(tree_path),
        method_anc: MethodAncestralCli::Marginal,
        dense: Some(dense),
        model_args: ModelArgs {
          model: GtrModelNameCli::JC69,
          ..ModelArgs::default()
        },
        impute_missing_data: true,
        output: OutputCoreArgs {
          output_tree_nwk: Some(dir.path().join("tree_out.nwk")),
          ..Default::default()
        },
        output_augur_node_data: Some(node_data_path.clone()),
        ..TreetimeAncestralArgsRaw::default()
      })
      .unwrap();

      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress).unwrap();
      json_read_str(std::fs::read_to_string(node_data_path).unwrap()).unwrap()
    }

    pub(super) fn reconstruct_json_with_translations() -> AugurNodeDataJsonAncestral {
      let dir = tempdir().unwrap();
      let tree_path = dir.path().join("tree.nwk");
      let fasta_path = dir.path().join("aln.fasta");
      let node_data_path = dir.path().join("augur-node-data.json");
      std::fs::write(&tree_path, "(A:0.1,B:0.1)root;").unwrap();
      std::fs::write(&fasta_path, ">A\nACGACG\n>B\nACGACG\n").unwrap();

      let translations_dir = dir.path().join("translations");
      std::fs::create_dir_all(&translations_dir).unwrap();
      std::fs::write(translations_dir.join("S.fasta"), ">A\nMKL\n>B\nMKL\n").unwrap();
      std::fs::write(translations_dir.join("M.fasta"), ">A\nWY\n>B\nWY\n").unwrap();
      let template = format!("{}/{{cds}}.fasta", translations_dir.display());

      let args = TreetimeAncestralArgs::try_from(TreetimeAncestralArgsRaw {
        alignment: AlignmentArgs {
          alignment: vec![fasta_path],
        },
        tree: Some(tree_path),
        method_anc: MethodAncestralCli::Marginal,
        dense: Some(false),
        model_args: ModelArgs {
          model: GtrModelNameCli::JC69,
          ..ModelArgs::default()
        },
        translations: Some(template),
        cdses: vec![o!("S"), o!("M")],
        output: OutputCoreArgs {
          output_tree_nwk: Some(dir.path().join("tree_out.nwk")),
          ..Default::default()
        },
        output_augur_node_data: Some(node_data_path.clone()),
        ..TreetimeAncestralArgsRaw::default()
      })
      .unwrap();

      run_ancestral_reconstruction(&args, &NoopCancel, &NoopProgress, &NoopProgress).unwrap();
      json_read_str(std::fs::read_to_string(node_data_path).unwrap()).unwrap()
    }

    pub(super) fn build_json_with_aa() -> AugurNodeDataJsonAncestral {
      let (graph, names, maps) = mutation_case();
      let name_to_key = node_name_to_key(&names, &graph);
      let mut aa_node_data = AaNodeData::default();

      aa_node_data.add_cds(
        "S",
        AaCdsNodeData {
          reference: "AC".to_owned(),
          root_sequence: "AC".to_owned(),
          node_mutations: btreemap! {
            name_to_key["A"] => vec![MutationEvent::Substitution(sub(b'C', 1, b'D'))],
            name_to_key["B"] => vec![],
            name_to_key["root"] => vec![],
          },
        },
      );

      let aa_annotations = btreemap! {
        o!("S") => AugurNodeDataJsonAnnotationEntry {
          start: Some(1),
          end: Some(6),
          strand: Some("+".to_owned()),
          entry_type: Some("CDS".to_owned()),
          segments: None,
          other: btreemap! {},
        },
      };

      let amino_acids = TreeAminoAcids {
        node_data: &aa_node_data,
        cdses: &aa_annotations,
      };
      build_json(&graph, &names, maps, &[false, false, false, false], Some(amino_acids))
    }

    pub(super) fn expected_json_with_aa() -> AugurNodeDataJsonAncestral {
      AugurNodeDataJsonAncestral {
        generated_by: Some(AugurNodeDataJsonGeneratedBy {
          program: "treetime".to_owned(),
          version: env!("CARGO_PKG_VERSION").to_owned(),
        }),
        metadata: AugurNodeDataJsonAncestralMeta {
          annotations: Some(AugurNodeDataJsonAnnotations {
            nuc: Some(AugurNodeDataJsonAnnotationEntry {
              start: Some(1),
              end: Some(4),
              strand: Some("+".to_owned()),
              entry_type: Some("source".to_owned()),
              segments: None,
              other: btreemap! {},
            }),
            other: btreemap! {
              o!("S") => AugurNodeDataJsonAnnotationEntry {
                start: Some(1),
                end: Some(6),
                strand: Some("+".to_owned()),
                entry_type: Some("CDS".to_owned()),
                segments: None,
                other: btreemap! {},
              },
            },
          }),
          reference: Some(btreemap! {
            o!("S") => o!("AC"),
            o!("nuc") => o!("ACGT"),
          }),
          mask: Some("0000".to_owned()),
          other: btreemap! {},
        },
        nodes: btreemap! {
          o!("A") => AugurNodeDataJsonAncestralNode {
            muts: vec![o!("T4A")],
            sequence: Some(o!("ACGA")),
            aa_muts: Some(btreemap! {
              o!("S") => vec![o!("C2D")],
            }),
            aa_sequences: None,
            other: btreemap! {},
          },
          o!("B") => AugurNodeDataJsonAncestralNode {
            muts: vec![],
            sequence: Some(o!("ACGT")),
            aa_muts: Some(btreemap! {
              o!("S") => vec![],
            }),
            aa_sequences: None,
            other: btreemap! {},
          },
          o!("root") => AugurNodeDataJsonAncestralNode {
            muts: vec![],
            sequence: Some(o!("ACGT")),
            aa_muts: Some(btreemap! {
              o!("S") => vec![],
            }),
            aa_sequences: Some(btreemap! {
              o!("S") => o!("AC"),
            }),
            other: btreemap! {},
          },
        },
      }
    }

    pub(super) fn expected_invariant_json() -> String {
      format!(
        r#"{{
  "generated_by": {{
    "program": "treetime",
    "version": "{version}"
  }},
  "nodes": {{
    "A": {{
      "muts": [],
      "sequence": "ACGT"
    }},
    "B": {{
      "muts": [],
      "sequence": "ACGT"
    }},
    "root": {{
      "muts": [],
      "sequence": "ACGT"
    }}
  }},
  "annotations": {{
    "nuc": {{
      "start": 1,
      "end": 4,
      "strand": "+",
      "type": "source"
    }}
  }},
  "reference": {{
    "nuc": "ACGT"
  }},
  "mask": "0000"
}}"#,
        version = env!("CARGO_PKG_VERSION")
      )
    }
  }
}
