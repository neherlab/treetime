#[cfg(test)]
mod tests {
  use pretty_assertions::assert_eq;
  use treetime_utils::io::json::json_read_str;
  use util_augur_node_data_json::AugurNodeDataJsonAncestral;

  // --- Full output, hand-built Fitch partition with one mutation ---

  #[test]
  fn test_augur_node_data_ancestral_full_output() {
    let (graph, partition) = helpers::mutation_case();
    let actual = helpers::write_json(&graph, &partition, &[false, false, false, false]);

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
    let (graph, partition) = helpers::mutation_case();
    let json_str = helpers::write_json(&graph, &partition, &[false, false, false, false]);

    let original: serde_json::Value = serde_json::from_str(&json_str).unwrap();
    let typed: AugurNodeDataJsonAncestral = json_read_str(&json_str).unwrap();
    let roundtripped: serde_json::Value = serde_json::to_value(&typed).unwrap();

    assert_eq!(original, roundtripped);
  }

  // --- Mask removes mutations at masked positions and masks output sequences,
  //     while the reference (root) sequence stays unmasked ---

  #[test]
  fn test_augur_node_data_ancestral_mask_filters_mutations() {
    let (graph, partition) = helpers::mutation_case();
    let actual = helpers::write_json(&graph, &partition, &[false, false, false, true]);

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
    let (graph, partition) = helpers::mutation_case();
    let json_str = helpers::write_json(&graph, &partition, &[false, false, false, false]);
    let data: AugurNodeDataJsonAncestral = json_read_str(&json_str).unwrap();

    // Root has no parent edge: empty mutations. A non-root node carries the
    // single mutation, confirming the empty-root assertion is not vacuous.
    assert_eq!(Vec::<String>::new(), data.nodes["root"].muts);
    assert_eq!(vec!["T4A".to_owned()], data.nodes["A"].muts);
  }

  // --- End-to-end through each reconstruction path with identical leaf
  //     sequences (deterministic: no mutations, root reconstructs to the
  //     shared sequence). Exercises the partition trait impls and the
  //     run.rs wiring for parsimony, marginal sparse, and marginal dense. ---

  #[test]
  fn test_augur_node_data_ancestral_parsimony_end_to_end() {
    use crate::ancestral::params::MethodAncestral;
    use crate::gtr::get_gtr::GtrModelName;
    let actual = helpers::reconstruct_json(MethodAncestral::Parsimony, None, GtrModelName::Infer);
    assert_eq!(helpers::expected_invariant_json(), actual.trim());
  }

  #[test]
  fn test_augur_node_data_ancestral_marginal_sparse_end_to_end() {
    use crate::ancestral::params::MethodAncestral;
    use crate::gtr::get_gtr::GtrModelName;
    let actual = helpers::reconstruct_json(MethodAncestral::Marginal, Some(false), GtrModelName::JC69);
    assert_eq!(helpers::expected_invariant_json(), actual.trim());
  }

  #[test]
  fn test_augur_node_data_ancestral_marginal_dense_end_to_end() {
    use crate::ancestral::params::MethodAncestral;
    use crate::gtr::get_gtr::GtrModelName;
    let actual = helpers::reconstruct_json(MethodAncestral::Marginal, Some(true), GtrModelName::JC69);
    assert_eq!(helpers::expected_invariant_json(), actual.trim());
  }

  mod helpers {
    use crate::alphabet::alphabet::Alphabet;
    use crate::ancestral::params::MethodAncestral;
    use crate::commands::ancestral::args::TreetimeAncestralArgs;
    use crate::commands::ancestral::augur_node_data::write_augur_node_data_json;
    use crate::commands::ancestral::run::run_ancestral_reconstruction;
    use crate::gtr::get_gtr::GtrModelName;
    use crate::partition::fitch::PartitionFitch;
    use crate::partition::sparse::{SparseEdgePartition, SparseNodePartition};
    use crate::payload::ancestral::GraphAncestral;
    use crate::progress::NoopProgress;
    use crate::seq::mutation::Sub;
    use maplit::btreemap;
    use std::collections::BTreeMap;
    use tempfile::{NamedTempFile, tempdir};
    use treetime_graph::node::Named;
    use treetime_io::nwk::nwk_read_str;
    use treetime_primitives::{AsciiChar, Seq};
    use treetime_utils::o;

    pub fn sub(reff: u8, pos: usize, qry: u8) -> Sub {
      Sub::new(
        AsciiChar::from_byte_unchecked(reff),
        pos,
        AsciiChar::from_byte_unchecked(qry),
      )
      .unwrap()
    }

    /// Two-leaf tree where leaf A differs from the root at one position.
    /// Root sequence ACGT, A is ACGA, with the substitution T4A on edge root->A.
    pub fn mutation_case() -> (GraphAncestral, PartitionFitch) {
      let graph: GraphAncestral = nwk_read_str("(A:0.1,B:0.1)root;").unwrap();
      let seqs = btreemap! { o!("A") => o!("ACGA"), o!("B") => o!("ACGT"), o!("root") => o!("ACGT") };
      let edge_subs = btreemap! { o!("A") => vec![sub(b'T', 3, b'A')] };
      let partition = build_fitch_partition(&graph, &seqs, &edge_subs, 4);
      (graph, partition)
    }

    pub fn build_fitch_partition(
      graph: &GraphAncestral,
      seqs: &BTreeMap<String, String>,
      edge_subs_by_child: &BTreeMap<String, Vec<Sub>>,
      length: usize,
    ) -> PartitionFitch {
      let alphabet = Alphabet::default();

      let mut key_to_name = BTreeMap::new();
      let mut nodes = BTreeMap::new();
      for node in graph.get_nodes() {
        let node_guard = node.read_arc();
        let key = node_guard.key();
        let payload = node_guard.payload().read_arc();
        let name = payload.name().unwrap().as_ref().to_owned();
        let seq = Seq::try_from_str(&seqs[&name]).unwrap();
        nodes.insert(key, SparseNodePartition::new(&seq, &alphabet).unwrap());
        key_to_name.insert(key, name);
      }

      let mut edges = BTreeMap::new();
      for edge in graph.get_edges() {
        let edge_guard = edge.read_arc();
        let edge_key = edge_guard.key();
        let child_name = &key_to_name[&edge_guard.target()];
        let subs = edge_subs_by_child.get(child_name).cloned().unwrap_or_default();
        edges.insert(edge_key, SparseEdgePartition::with_fitch_subs(subs));
      }

      PartitionFitch {
        index: 0,
        alphabet,
        length,
        nodes,
        edges,
      }
    }

    pub fn write_json(graph: &GraphAncestral, partition: &PartitionFitch, mask: &[bool]) -> String {
      let tmp = NamedTempFile::new().unwrap();
      write_augur_node_data_json(graph, partition, mask, tmp.path()).unwrap();
      std::fs::read_to_string(tmp.path()).unwrap()
    }

    pub fn reconstruct_json(method: MethodAncestral, dense: Option<bool>, model: GtrModelName) -> String {
      let dir = tempdir().unwrap();
      let tree_path = dir.path().join("tree.nwk");
      let fasta_path = dir.path().join("aln.fasta");
      std::fs::write(&tree_path, "(A:0.1,B:0.1)root;").unwrap();
      std::fs::write(&fasta_path, ">A\nACGT\n>B\nACGT\n").unwrap();

      let args = TreetimeAncestralArgs {
        input_fastas: vec![fasta_path],
        tree: tree_path,
        method_anc: method,
        dense,
        model_name: model,
        outdir: dir.path().to_path_buf(),
        ..TreetimeAncestralArgs::default()
      };

      run_ancestral_reconstruction(&args, &NoopProgress).unwrap();
      std::fs::read_to_string(dir.path().join("ancestral.augur-node-data.json")).unwrap()
    }

    /// Expected JSON when all leaves share the same sequence ACGT: every node
    /// reconstructs to ACGT, no mutations anywhere, empty mask.
    pub fn expected_invariant_json() -> String {
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
