#[cfg(test)]
mod tests {
  use crate::commands::clock::clock_graph::GraphClock;
  use crate::graph::__tests__::graph::tests::{TestEdge, TestNode};
  use crate::o;
  use crate::seq::div::{OnlyLeaves, compute_divs};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_calculate_divs() -> Result<(), Report> {
    let graph: Graph<TestNode, TestEdge, ()> = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = compute_divs(&graph, OnlyLeaves(false));
    let expected = btreemap! {
      o!("A") => 0.20000000298023224,
      o!("AB") => 0.10000000149011612,
      o!("B") => 0.30000000447034836,
      o!("C") => 0.2500000037252903,
      o!("CD") => 0.05000000074505806,
      o!("D") => 0.16999999806284904,
      o!("root") => 0.0,
    };
    assert_eq!(&expected, &actual);
    Ok(())
  }

  #[test]
  fn test_calculate_divs_only_leaves() -> Result<(), Report> {
    let graph: GraphClock = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let actual = compute_divs(&graph, OnlyLeaves(true));
    let expected = btreemap! {
      o!("A") => 0.20000000298023224,
      o!("B") => 0.30000000447034836,
      o!("C") => 0.2500000037252903,
      o!("D") => 0.16999999806284904,
    };
    assert_eq!(&expected, &actual);
    Ok(())
  }
}
