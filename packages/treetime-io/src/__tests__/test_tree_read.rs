#[cfg(test)]
mod tests {
  use crate::nwk::{NewickDialect, NwkNodeComments, NwkWriteOptions, TREE_DIALECT_DEFAULT, nwk_write_str};
  use crate::tree::tree_read;
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::tree_view::TreeView;

  #[rustfmt::skip]
  #[rstest]
  #[case::newick(            "(A:0.1,B:0.2)root;")]
  #[case::nexus(             "#NEXUS\nBegin Trees;\n  Tree t1 = (A:0.1,B:0.2)root;\nEnd;\n")]
  #[case::nexus_translate(   "#NEXUS\nBegin Trees;\n  Translate 1 A, 2 B;\n  Tree t1 = [&R] (1:0.1,2:0.2)root;\nEnd;\n")]
  #[case::nexus_lowercase(   "\u{feff}  #nexus\nbegin trees;\n  tree t1 = (A:0.1,B:0.2)root;\nend;\n")]
  #[trace]
  fn test_tree_read_newick_and_nexus(#[case] input: &str) -> Result<(), Report> {
    let parse = tree_read(input.as_bytes(), TREE_DIALECT_DEFAULT)?;

    let written = nwk_write_str(
      &TreeView::new(&parse.graph)?,
      &parse.names(),
      &parse.branch_lengths,
      &NwkWriteOptions::default(),
      &NwkNodeComments::new(),
    )?;
    assert_eq!("(A:0.1,B:0.2)root;", written);
    Ok(())
  }

  #[test]
  fn test_tree_read_nexus_with_two_trees_is_error() {
    let input = indoc! {"
      #NEXUS
      Begin Trees;
        Tree t1 = (A,B);
        Tree t2 = (A,B);
      End;
    "};

    let actual = format!("{:#}", tree_read(input.as_bytes(), TREE_DIALECT_DEFAULT).unwrap_err());

    assert_eq!(
      "The Nexus file contains 2 trees, but TreeTime reads exactly one tree",
      actual
    );
  }

  #[test]
  fn test_tree_read_nexus_without_trees_is_error() {
    let actual = format!(
      "{:#}",
      tree_read(b"#NEXUS\nBegin Taxa;\nEnd;\n".as_slice(), TREE_DIALECT_DEFAULT).unwrap_err()
    );

    assert_eq!("The Nexus file contains no tree", actual);
  }

  #[test]
  fn test_tree_read_hash_in_names() -> Result<(), Report> {
    let parse = tree_read(
      "(A#1:0.1,'B#2':0.2,EPI_ISL#402124:0.3)root;".as_bytes(),
      TREE_DIALECT_DEFAULT,
    )?;

    let written = nwk_write_str(
      &TreeView::new(&parse.graph)?,
      &parse.names(),
      &parse.branch_lengths,
      &NwkWriteOptions::default(),
      &NwkNodeComments::new(),
    )?;
    assert_eq!("('A#1':0.1,'B#2':0.2,'EPI_ISL#402124':0.3)root;", written);
    Ok(())
  }

  #[test]
  fn test_tree_read_infinite_length_cannot_be_written() -> Result<(), Report> {
    let parse = tree_read(b"(A:1e400,B:1)root;".as_slice(), TREE_DIALECT_DEFAULT)?;

    let actual = nwk_write_str(
      &TreeView::new(&parse.graph)?,
      &parse.names(),
      &parse.branch_lengths,
      &NwkWriteOptions::default(),
      &NwkNodeComments::new(),
    );

    assert_eq!(
      "When writing Newick: When writing the branch above node 1 ('A'): Newick cannot represent the number inf",
      format!("{:#}", actual.unwrap_err())
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beast_keeps_hash_in_name(  NewickDialect::BEAST,   "('EPI_ISL#402124':0.3,B:0.2)root;")]
  #[case::enewick_reads_hybrid_tag(  NewickDialect::ENEWICK, "(EPI_ISL:0.3,B:0.2)root;")]
  #[trace]
  fn test_tree_read_hash_tag_by_dialect(#[case] dialect: NewickDialect, #[case] expected: &str) -> Result<(), Report> {
    let parse = tree_read(b"(EPI_ISL#402124:0.3,B:0.2)root;".as_slice(), dialect)?;

    let written = nwk_write_str(
      &TreeView::new(&parse.graph)?,
      &parse.names(),
      &parse.branch_lengths,
      &NwkWriteOptions::default(),
      &NwkNodeComments::new(),
    )?;
    assert_eq!(expected, written);
    Ok(())
  }

  #[test]
  fn test_tree_read_network_is_error() {
    let input = "((A[&segments={0,1}]:1,(B[&segments={0,1}]:1)#H0[&segments={0}]:0.5)[&segments={0,1}]:1,(#H0[&segments={1}]:0.7,C[&segments={0,1}]:1)[&segments={0,1}]:1);";

    let actual = format!(
      "{:#}",
      tree_read(input.as_bytes(), NewickDialect::ENEWICK_BEAST).unwrap_err()
    );

    assert_eq!(
      "When reading Newick: The tree contains the network node '#H0' with 2 parents; TreeTime reads trees only",
      actual
    );
  }
}
