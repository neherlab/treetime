#[cfg(test)]
mod tests {
  use crate::nex::{NexWriteOptions, nex_write_str};
  use crate::nwk::{EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, nwk_read_str};
  use eyre::Report;
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::{GraphEdge, HasBranchLength};
  use treetime_graph::node::{GraphNode, Named};

  #[derive(Clone, Debug, Default)]
  struct TestNode {
    name: Option<String>,
  }

  impl GraphNode for TestNode {}

  impl Named for TestNode {
    fn name(&self) -> Option<impl AsRef<str>> {
      self.name.as_deref()
    }

    fn set_name(&mut self, name: Option<impl AsRef<str>>) {
      self.name = name.map(|n| n.as_ref().to_owned());
    }
  }

  impl NodeFromNwk for TestNode {
    fn from_nwk(
      name: Option<impl AsRef<str>>,
      _confidence: Option<f64>,
      _: &BTreeMap<String, String>,
    ) -> Result<Self, Report> {
      Ok(Self {
        name: name.map(|n| n.as_ref().to_owned()),
      })
    }
  }

  impl NodeToNwk for TestNode {
    fn nwk_name(&self) -> Option<impl AsRef<str>> {
      self.name.as_deref()
    }
  }

  #[derive(Clone, Debug, Default)]
  struct TestEdge(Option<f64>);

  impl GraphEdge for TestEdge {}

  impl HasBranchLength for TestEdge {
    fn branch_length(&self) -> Option<f64> {
      self.0
    }

    fn set_branch_length(&mut self, weight: Option<f64>) {
      self.0 = weight;
    }
  }

  impl EdgeFromNwk for TestEdge {
    fn from_nwk(weight: Option<f64>) -> Result<Self, Report> {
      Ok(Self(weight))
    }
  }

  impl EdgeToNwk for TestEdge {
    fn nwk_weight(&self) -> Option<f64> {
      self.0
    }
  }

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
    let graph = nwk_read_str::<TestNode, TestEdge, ()>(nwk)?;
    let actual = nex_write_str(&graph, &NexWriteOptions::default())?;
    assert_eq!(expected, actual);
    Ok(())
  }
}
