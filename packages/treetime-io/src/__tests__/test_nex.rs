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
}
