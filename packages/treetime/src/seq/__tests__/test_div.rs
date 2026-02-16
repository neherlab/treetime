#[cfg(test)]
mod tests {
  use crate::commands::clock::clock_graph::GraphClock;
  use crate::graph::__tests__::graph::tests::{TestEdge, TestNode};
  use crate::o;
  use crate::seq::div::{OnlyLeaves, compute_divs};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_calculate_divs() -> Result<(), Report> {
    let graph: Graph<TestNode, TestEdge, ()> = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = compute_divs(&graph, OnlyLeaves(false));

    // Expected values are exact sums of branch lengths from root to each node
    // Tree: ((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01
    let expected = btreemap! {
      o!("A") => 0.2,      // root->AB (0.1) + AB->A (0.1)
      o!("AB") => 0.1,     // root->AB (0.1)
      o!("B") => 0.3,      // root->AB (0.1) + AB->B (0.2)
      o!("C") => 0.25,     // root->CD (0.05) + CD->C (0.2)
      o!("CD") => 0.05,    // root->CD (0.05)
      o!("D") => 0.17,     // root->CD (0.05) + CD->D (0.12)
      o!("root") => 0.0,
    };

    assert_eq!(expected.len(), actual.len());
    for (name, expected_div) in &expected {
      assert_abs_diff_eq!(actual[name], *expected_div, epsilon = 1e-6);
    }

    Ok(())
  }

  #[test]
  fn test_calculate_divs_only_leaves() -> Result<(), Report> {
    let graph: GraphClock = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let actual = compute_divs(&graph, OnlyLeaves(true));

    let expected = btreemap! {
      o!("A") => 0.2,
      o!("B") => 0.3,
      o!("C") => 0.25,
      o!("D") => 0.17,
    };

    assert_eq!(expected.len(), actual.len());
    for (name, expected_div) in &expected {
      assert_abs_diff_eq!(actual[name], *expected_div, epsilon = 1e-6);
    }

    Ok(())
  }
}
