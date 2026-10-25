#[cfg(test)]
mod tests {
  use crate::__tests__::test_support::tests::RecordingLog;
  use crate::commands::mugration::run::pair_traits;
  use itertools::Itertools;
  use maplit::btreeset;
  use pretty_assertions::assert_eq;
  use std::path::Path;
  use treetime::mugration::trait_pairing::TraitPairing;
  use treetime::progress::{RunWarning, RunWarningKind};
  use treetime_graph::keymap;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::o;

  #[test]
  fn test_pair_traits_keeps_the_first_row_of_a_name_for_every_leaf_with_that_name() {
    let tree = nwk_read(b"((A:0.1,A:0.1)X:0.1,B:0.1)root;".as_slice()).unwrap();
    let names = tree.names.clone();
    let rows = vec![
      (o!("A"), o!("usa")),
      (o!("B"), o!("uk")),
      (o!("A"), o!("france")),
      (o!("Z"), o!("germany")),
    ];
    let log = RecordingLog::default();

    let pairing = pair_traits(rows, &tree.graph, &names, Path::new("metadata.tsv"), &log);

    let leaf_keys = |name: &str| {
      tree
        .graph
        .get_leaves()
        .map(|leaf| leaf.key())
        .filter(|key| names.pairing_name(*key) == name)
        .sorted()
        .collect_vec()
    };
    let [a1, a2] = leaf_keys("A")[..] else {
      panic!("the tree has two leaves named A")
    };
    assert_eq!(
      (
        TraitPairing {
          traits: keymap! { a1 => o!("usa"), a2 => o!("usa"), leaf_keys("B")[0] => o!("uk") },
          observed_values: btreeset! { o!("germany"), o!("uk"), o!("usa") },
          duplicate_names: vec![o!("A")],
          unmatched: vec![(o!("Z"), o!("germany"))],
        },
        vec![RunWarning {
          kind: RunWarningKind::DuplicateMetadataNames,
          message: o!(
            "The metadata 'metadata.tsv' has more than one row named A. TreeTime uses the first row of each name."
          ),
          names: vec![o!("A")],
        }],
      ),
      (pairing, log.into_warnings())
    );
  }
}
