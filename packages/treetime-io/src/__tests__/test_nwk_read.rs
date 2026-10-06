#[cfg(test)]
mod tests {
  use crate::nwk::{NwkParse, NwkWriteOptions, nwk_read, nwk_write_str};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::tree_view::TreeView;

  #[rustfmt::skip]
  #[rstest]
  #[case::unique_names(          "(A:0.1,B:0.1,C:0.1)root;",          vec![])]
  #[case::repeated_leaf(         "(A:0.1,A:0.1,B:0.1)root;",          vec!["A"])]
  #[case::leaf_and_internal_node("((A:0.1,B:0.1)A:0.1,B:0.1)root;",   vec!["A", "B"])]
  #[case::unnamed_nodes(         "((A:0.1,B:0.1):0.1,C:0.1);",        vec![])]
  #[trace]
  fn test_nwk_read_lists_repeated_node_names(#[case] input: &str, #[case] expected: Vec<&str>) -> Result<(), Report> {
    let parse = nwk_read(input.as_bytes())?;

    assert_eq!(expected, parse.duplicate_names);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::integer_support( "((A:0.1,B:0.2)95:0.3,C:0.4);",    "((A:0.1,B:0.2)NODE_0000001:0.3,C:0.4)NODE_0000000;")]
  #[case::fraction_support("((A:0.1,B:0.2)0.959:0.3,C:0.4);", "((A:0.1,B:0.2)NODE_0000001:0.3,C:0.4)NODE_0000000;")]
  #[case::name(            "((A:0.1,B:0.2)inner:0.3,C:0.4);",  "((A:0.1,B:0.2)inner:0.3,C:0.4)NODE_0000000;")]
  #[trace]
  fn test_nwk_read_internal_label(#[case] input: &str, #[case] expected: &str) -> Result<(), Report> {
    let parse = nwk_read(input.as_bytes())?;
    let names = parse.names();
    let NwkParse {
      graph, branch_lengths, ..
    } = parse;

    let actual = nwk_write_str(
      &TreeView::new(&graph)?,
      &names,
      &branch_lengths,
      &NwkWriteOptions::default(),
      &btreemap! {},
    )?;

    assert_eq!(expected, actual);
    Ok(())
  }
}
