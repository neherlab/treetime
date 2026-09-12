#[cfg(test)]
mod tests {
  use crate::nex::{NexWriteOptions, nex_write_str};
  use crate::nwk::{NwkParse, nwk_read_str};
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

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
  #[trace]
  fn test_nex_exact_output(#[case] nwk: &str, #[case] expected: &str) -> Result<(), Report> {
    let NwkParse {
      graph,
      names,
      branch_lengths,
      ..
    } = nwk_read_str::<()>(nwk)?;
    let actual = nex_write_str(&graph, &names, &branch_lengths, &NexWriteOptions::default())?;
    assert_eq!(expected, actual);
    Ok(())
  }
}
