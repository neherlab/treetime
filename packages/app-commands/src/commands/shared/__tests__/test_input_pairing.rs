#[cfg(test)]
mod tests {
  use crate::commands::shared::alignment::{PairedAlignment, pair_alignment};
  use crate::commands::shared::dates_input::read_input_dates;
  use crate::commands::shared::metadata::MetadataIdArgs;
  use crate::commands::shared::topology_order_args::target_positions;
  use crate::commands::shared::tree_input::read_input_tree;
  use crate::runs::warnings::WarningCollector;
  use helpers::{leaf_keys, warning};
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::fs;
  use std::path::PathBuf;
  use tempfile::tempdir;
  use treetime::progress::{NoopProgress, RunWarningKind};
  use treetime_io::dates_csv::DateConstraint;
  use treetime_io::nwk::{TREE_DIALECT_DEFAULT, nwk_read};
  use treetime_primitives::Seq;
  use treetime_utils::{o, vec_of_owned};

  const TREE: &str = "((A:0.1,A:0.1)X:0.1,B:0.1)root;";

  #[test]
  fn test_input_pairing_tree_with_a_repeated_name_raises_one_warning() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("tree.nwk");
    fs::write(&path, format!("{TREE}\n")).unwrap();
    let log = WarningCollector::new(&NoopProgress);

    read_input_tree(&path, TREE_DIALECT_DEFAULT, &log).unwrap();

    assert_eq!(
      vec![warning(
        RunWarningKind::DuplicateNodeNames,
        &format!(
          "The tree '{}' gives the same name to more than one node: A. Nodes with the same name receive the same data from the other inputs, and augur node data keeps one entry per name.",
          path.display()
        ),
      )],
      log.into_warnings()
    );
  }

  #[test]
  fn test_input_pairing_alignment_keeps_the_first_record_of_a_name_for_every_leaf_with_that_name() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let names = tree.names();
    let records = vec![
      helpers::record("A", Some("first"), "ACGT"),
      helpers::record("B", None, "GGGG"),
      helpers::record("A", Some("second"), "TTTT"),
      helpers::record("Z", None, "CCCC"),
    ];
    let log = WarningCollector::new(&NoopProgress);

    let PairedAlignment { sequences, descs } =
      pair_alignment(records, &[PathBuf::from("aln.fasta")], &tree.graph, &names, &log);

    let [a1, a2] = leaf_keys(&tree, "A")[..] else {
      panic!("the tree has two leaves named A")
    };
    let b = leaf_keys(&tree, "B")[0];
    let leaf_seqs = [a1, a2, b].map(|key| sequences.nodes[&key].seq.as_ref().map(ToString::to_string));
    assert_eq!(
      (
        [Some(o!("ACGT")), Some(o!("ACGT")), Some(o!("GGGG"))],
        btreemap! { a1 => Some(o!("first")), a2 => Some(o!("first")), b => None },
        vec![(o!("Z"), Seq::try_from_str("CCCC").unwrap())],
        vec![warning(
          RunWarningKind::DuplicateSequenceNames,
          "The alignment 'aln.fasta' has more than one sequence named A. TreeTime uses the first sequence of each name.",
        )],
      ),
      (leaf_seqs, descs, sequences.unmatched, log.into_warnings())
    );
  }

  #[test]
  fn test_input_pairing_dates_keep_the_first_row_of_a_name_for_every_node_with_that_name() {
    let dir = tempdir().unwrap();
    let path = dir.path().join("metadata.tsv");
    fs::write(&path, "strain\tdate\nA\t2001\nB\t2002\nA\t2005\n").unwrap();
    let tree = nwk_read(TREE.as_bytes()).unwrap();
    let names = tree.names();
    let log = WarningCollector::new(&NoopProgress);

    let actual = read_input_dates(&path, &MetadataIdArgs::default(), None, &tree.graph, &names, &log).unwrap();

    let [a1, a2] = leaf_keys(&tree, "A")[..] else {
      panic!("the tree has two leaves named A")
    };
    let b = leaf_keys(&tree, "B")[0];
    let date = |value: &str| DateConstraint {
      raw: value.to_owned(),
      ..DateConstraint::exact(value.parse().unwrap())
    };
    assert_eq!(
      (
        btreemap! { a1 => date("2001"), a2 => date("2001"), b => date("2002") },
        vec![warning(
          RunWarningKind::DuplicateMetadataNames,
          &format!(
            "The metadata '{}' has more than one row named A. TreeTime uses the first row of each name.",
            path.display()
          ),
        )],
      ),
      (actual, log.into_warnings())
    );
  }

  #[test]
  fn test_input_pairing_target_order_gives_leaves_that_share_a_name_its_first_position() {
    let tree = nwk_read(TREE.as_bytes()).unwrap();

    let actual = target_positions(vec_of_owned!["B", "A", "Z", "A"], &tree.graph, &tree.names());

    let [a1, a2] = leaf_keys(&tree, "A")[..] else {
      panic!("the tree has two leaves named A")
    };
    assert_eq!(btreemap! { leaf_keys(&tree, "B")[0] => 0, a1 => 1, a2 => 1 }, actual);
  }

  mod helpers {
    use itertools::Itertools;
    use treetime::progress::{RunWarning, RunWarningKind};
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::fasta::FastaRecord;
    use treetime_io::nwk::NwkParse;
    use treetime_primitives::Seq;
    use treetime_utils::o;

    pub(super) fn leaf_keys(tree: &NwkParse, name: &str) -> Vec<GraphNodeKey> {
      let names = tree.names();
      tree
        .graph
        .get_leaves()
        .map(|leaf| leaf.key())
        .filter(|key| names[key].as_deref() == Some(name))
        .sorted()
        .collect()
    }

    pub(super) fn record(name: &str, desc: Option<&str>, seq: &str) -> FastaRecord {
      FastaRecord {
        seq_name: name.to_owned(),
        desc: desc.map(str::to_owned),
        seq: Seq::try_from_str(seq).unwrap(),
      }
    }

    pub(super) fn warning(kind: RunWarningKind, message: &str) -> RunWarning {
      RunWarning {
        kind,
        message: message.to_owned(),
        names: vec![o!("A")],
      }
    }
  }
}
