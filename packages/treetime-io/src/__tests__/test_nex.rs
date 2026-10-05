#[cfg(test)]
mod tests {
  use crate::nex::nex_write;
  use crate::nwk::{NwkParse, NwkWriteOptions, nwk_read};
  use eyre::Report;
  use indoc::indoc;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::tree_view::TreeView;

  #[rustfmt::skip]
  #[rstest]
  #[case::two_leaves(
    "(A:0.1,B:0.2)root;",
    indoc! {r#"
      #NEXUS
      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels A B;
      End;
      Begin Trees;
        Tree tree1=(A:0.1,B:0.2)root;
      End;
    "#},
  )]
  #[case::nested_tree(
    "((C:0.3,D:0.4)E:0.5,F:0.1)G;",
    indoc! {r#"
      #NEXUS
      Begin Taxa;
        Dimensions NTax=3;
        TaxLabels C D F;
      End;
      Begin Trees;
        Tree tree1=((C:0.3,D:0.4)E:0.5,F:0.1)G;
      End;
    "#},
  )]
  #[case::single_leaf(
    "(A:0.1)root;",
    indoc! {r#"
      #NEXUS
      Begin Taxa;
        Dimensions NTax=1;
        TaxLabels A;
      End;
      Begin Trees;
        Tree tree1=(A:0.1)root;
      End;
    "#},
  )]
  #[case::quoted_labels(
    "('A B':0.1,'C''D':0.2)root;",
    indoc! {r#"
      #NEXUS
      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels 'A B' 'C''D';
      End;
      Begin Trees;
        Tree tree1=('A B':0.1,'C''D':0.2)root;
      End;
    "#},
  )]
  #[trace]
  fn test_nex_exact_output(#[case] nwk: &str, #[case] expected: &str) -> Result<(), Report> {
    let parse = nwk_read(nwk.as_bytes())?;
    let names = parse.names();
    let NwkParse {
      graph, branch_lengths, ..
    } = parse;
    let mut actual = Vec::new();
    nex_write(
      &mut actual,
      &TreeView::new(&graph)?,
      &names,
      &branch_lengths,
      &NwkWriteOptions::default(),
      &btreemap! {},
    )?;
    assert_eq!(expected, String::from_utf8(actual)?);
    Ok(())
  }

  #[test]
  fn test_nex_taxon_count_matches_written_labels() -> Result<(), Report> {
    let parse = nwk_read(b"(A:0.1,B:0.2)root;".as_slice())?;
    let mut names = parse.names();
    let b_key = names
      .iter()
      .find_map(|(key, name)| (name.as_deref() == Some("B")).then_some(*key))
      .expect("fixture has leaf B");
    names.insert(b_key, None);
    let NwkParse {
      graph, branch_lengths, ..
    } = parse;
    let mut actual = Vec::new();

    nex_write(
      &mut actual,
      &TreeView::new(&graph)?,
      &names,
      &branch_lengths,
      &NwkWriteOptions::default(),
      &btreemap! {},
    )?;

    let expected = indoc! {r#"
      #NEXUS
      Begin Taxa;
        Dimensions NTax=1;
        TaxLabels A;
      End;
      Begin Trees;
        Tree tree1=(A:0.1,:0.2)root;
      End;
    "#};
    assert_eq!(expected, String::from_utf8(actual)?);
    Ok(())
  }
}
