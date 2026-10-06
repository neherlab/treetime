#[cfg(test)]
mod tests {
  use crate::commands::mugration::run::{MugrationTraits, pair_traits};
  use crate::runs::warnings::WarningCollector;
  use itertools::Itertools;
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use std::path::Path;
  use treetime::progress::{NoopProgress, RunWarning, RunWarningKind};
  use treetime_io::nwk::nwk_read;
  use treetime_utils::o;

  #[test]
  fn test_pair_traits_keeps_the_first_row_of_a_name_for_every_leaf_with_that_name() {
    let tree = nwk_read(b"((A:0.1,A:0.1)X:0.1,B:0.1)root;".as_slice()).unwrap();
    let names = tree.names();
    let rows = vec![
      (o!("A"), o!("usa")),
      (o!("B"), o!("uk")),
      (o!("A"), o!("france")),
      (o!("Z"), o!("germany")),
    ];
    let log = WarningCollector::new(&NoopProgress);

    let MugrationTraits {
      traits,
      observed_values,
    } = pair_traits(rows, &tree.graph, &names, Path::new("metadata.tsv"), &log);

    let leaf_keys = |name: &str| {
      tree
        .graph
        .get_leaves()
        .map(|leaf| leaf.key())
        .filter(|key| names[key].as_deref() == Some(name))
        .sorted()
        .collect_vec()
    };
    let [a1, a2] = leaf_keys("A")[..] else {
      panic!("the tree has two leaves named A")
    };
    assert_eq!(
      (
        btreemap! { a1 => o!("usa"), a2 => o!("usa"), leaf_keys("B")[0] => o!("uk") },
        btreeset! { o!("germany"), o!("uk"), o!("usa") },
        vec![RunWarning {
          kind: RunWarningKind::DuplicateMetadataNames,
          message: o!(
            "The metadata 'metadata.tsv' has more than one row named A. TreeTime uses the first row of each name."
          ),
          names: vec![o!("A")],
        }],
      ),
      (traits, observed_values, log.into_warnings())
    );
  }
}
