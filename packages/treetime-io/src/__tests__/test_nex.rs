#[cfg(test)]
mod tests {
  use crate::nex::nex_write;
  use crate::nwk::{NwkParse, NwkWriteOptions, nwk_read};
  use eyre::Report;
  use indoc::indoc;
  use treetime_graph::keymap;

  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::schedule::GraphSchedule;

  #[rustfmt::skip]
  #[rstest]
  #[case::two_leaves(
    "(A:0.1,B:0.2)root;",
    indoc! {r#"
      #NEXUS

      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels
          A
          B
        ;
      End;

      Begin Trees;
        Tree tree1 = (A:0.1,B:0.2)root;
      End;
    "#},
  )]
  #[case::nested_tree(
    "((C:0.3,D:0.4)E:0.5,F:0.1)G;",
    indoc! {r#"
      #NEXUS

      Begin Taxa;
        Dimensions NTax=3;
        TaxLabels
          C
          D
          F
        ;
      End;

      Begin Trees;
        Tree tree1 = ((C:0.3,D:0.4)E:0.5,F:0.1)G;
      End;
    "#},
  )]
  #[case::single_leaf(
    "(A:0.1)root;",
    indoc! {r#"
      #NEXUS

      Begin Taxa;
        Dimensions NTax=1;
        TaxLabels
          A
        ;
      End;

      Begin Trees;
        Tree tree1 = (A:0.1)root;
      End;
    "#},
  )]
  #[case::quoted_labels(
    "('A B':0.1,'C''D':0.2)root;",
    indoc! {r#"
      #NEXUS

      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels
          'A B'
          'C''D'
        ;
      End;

      Begin Trees;
        Tree tree1 = ('A B':0.1,'C''D':0.2)root;
      End;
    "#},
  )]
  #[case::nexus_punctuation(
    "('hCoV-19/USA/1':0.1,'hCoV-19/USA/2':0.2)root;",
    indoc! {r#"
      #NEXUS

      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels
          'hCoV-19/USA/1'
          'hCoV-19/USA/2'
        ;
      End;

      Begin Trees;
        Tree tree1 = (hCoV-19/USA/1:0.1,hCoV-19/USA/2:0.2)root;
      End;
    "#},
  )]
  #[trace]
  fn test_nex_exact_output(#[case] nwk: &str, #[case] expected: &str) -> Result<(), Report> {
    let parse = nwk_read(nwk.as_bytes())?;
    let names = parse.names.clone();
    let NwkParse {
      graph, branch_lengths, ..
    } = parse;
    let mut actual = Vec::new();
    nex_write(
      &mut actual,
      &GraphSchedule::new(&graph)?.into_tree()?,
      &names,
      &branch_lengths,
      &NwkWriteOptions::default(),
      &keymap! {},
    )?;
    assert_eq!(expected, String::from_utf8(actual)?);
    Ok(())
  }

  #[test]
  fn test_nex_writes_the_generated_name_of_an_unnamed_leaf() -> Result<(), Report> {
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read(b"(A:0.1,:0.2)root;".as_slice())?;
    let mut actual = Vec::new();

    nex_write(
      &mut actual,
      &GraphSchedule::new(&graph)?.into_tree()?,
      &names,
      &branch_lengths,
      &NwkWriteOptions::default(),
      &keymap! {},
    )?;

    let expected = indoc! {r#"
      #NEXUS

      Begin Taxa;
        Dimensions NTax=2;
        TaxLabels
          A
          NODE_0000000
        ;
      End;

      Begin Trees;
        Tree tree1 = (A:0.1,NODE_0000000:0.2)root;
      End;
    "#};
    assert_eq!(expected, String::from_utf8(actual)?);
    Ok(())
  }
}
