#[cfg(test)]
mod tests {
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::optimization::outliers::mark_outlier_leaves;
  use eyre::Report;
  use maplit::{btreemap, btreeset};
  use pretty_assertions::assert_eq;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_mark_outlier_leaves_adds_outliers_to_leaf_flags() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)AB:0.1,C:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");

    let outliers = btreeset! { key("B") };
    let leaf_bad_branches = btreemap! {
      key("A") => true,
      key("B") => false,
      key("C") => false,
    };

    let actual = mark_outlier_leaves(&graph, &outliers, &leaf_bad_branches);

    let expected = btreemap! {
      key("A") => true,
      key("B") => true,
      key("C") => false,
    };
    assert_eq!(expected, actual);

    Ok(())
  }
}
