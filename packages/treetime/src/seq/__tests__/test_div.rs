#[cfg(test)]
mod tests {
  use crate::commands::clock::clock_graph::GraphClock;
  use crate::graph::__tests__::graph::tests::{TestEdge, TestNode};
  use crate::o;
  use crate::seq::div::{OnlyLeaves, compute_divs};
  use approx::assert_abs_diff_eq;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use treetime_graph::graph::Graph;
  use treetime_io::nwk::nwk_read_str;

  // OnlyLeaves(false) - all nodes
  #[test]
  fn test_all_nodes() -> Result<(), Report> {
    let graph: Graph<TestNode, TestEdge, ()> = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;

    let actual = compute_divs(&graph, OnlyLeaves(false));

    let expected = btreemap! {
      o!("root") => 0.0,
      o!("AB") => 0.1,
      o!("CD") => 0.05,
      o!("A") => 0.2,
      o!("B") => 0.3,
      o!("C") => 0.25,
      o!("D") => 0.17,
    };

    assert_eq!(expected.len(), actual.len());
    for (name, expected_div) in &expected {
      assert_abs_diff_eq!(actual[name], *expected_div, epsilon = 1e-8);
    }

    Ok(())
  }

  // OnlyLeaves(true) - leaves only
  #[test]
  fn test_only_leaves() -> Result<(), Report> {
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
      assert_abs_diff_eq!(actual[name], *expected_div, epsilon = 1e-8);
    }

    Ok(())
  }

  // Unnamed internal nodes get auto-generated names
  #[test]
  fn test_unnamed_internals() -> Result<(), Report> {
    let graph: Graph<TestNode, TestEdge, ()> = nwk_read_str("((A:0.1,B:0.2):0.1,(C:0.2,D:0.12):0.05):0.01;")?;

    let actual = compute_divs(&graph, OnlyLeaves(true));

    let expected = btreemap! {
      o!("A") => 0.2,
      o!("B") => 0.3,
      o!("C") => 0.25,
      o!("D") => 0.17,
    };

    assert_eq!(expected.len(), actual.len());
    for (name, expected_div) in &expected {
      assert_abs_diff_eq!(actual[name], *expected_div, epsilon = 1e-8);
    }

    Ok(())
  }

  // Single node tree
  #[test]
  fn test_single_node() -> Result<(), Report> {
    let graph: Graph<TestNode, TestEdge, ()> = nwk_read_str("A:0.5;")?;

    let actual = compute_divs(&graph, OnlyLeaves(true));

    assert_eq!(1, actual.len());
    assert_abs_diff_eq!(actual["A"], 0.0, epsilon = 1e-9);

    Ok(())
  }

  // Linear chain (no branching)
  #[test]
  fn test_linear_chain() -> Result<(), Report> {
    let graph: Graph<TestNode, TestEdge, ()> = nwk_read_str("((A:0.1)B:0.2)C:0.3;")?;

    let actual = compute_divs(&graph, OnlyLeaves(true));

    assert_eq!(1, actual.len());
    assert_abs_diff_eq!(actual["A"], 0.3, epsilon = 1e-8);

    Ok(())
  }

  // Deep tree (20 levels)
  #[test]
  fn test_deep_tree() -> Result<(), Report> {
    let depth = 20;
    let branch_len = 0.05;

    let mut nwk = format!("A:{branch_len}");
    for _ in 1..depth {
      nwk = format!("({nwk}):{branch_len}");
    }
    nwk.push(';');

    let graph: Graph<TestNode, TestEdge, ()> = nwk_read_str(&nwk)?;
    let actual = compute_divs(&graph, OnlyLeaves(true));

    assert_eq!(1, actual.len());
    let expected = (depth - 1) as f64 * branch_len;
    // 19 successive f64 additions of inexact 0.05; measured error 1.42e-8.
    assert_abs_diff_eq!(actual["A"], expected, epsilon = 1e-7);

    Ok(())
  }

  // Zero branch lengths
  #[test]
  fn test_zero_branch_lengths() -> Result<(), Report> {
    let graph: Graph<TestNode, TestEdge, ()> = nwk_read_str("((A:0.0,B:0.1):0.0,(C:0.2,D:0.0):0.1):0.0;")?;

    let actual = compute_divs(&graph, OnlyLeaves(true));

    let expected = btreemap! {
      o!("A") => 0.0,
      o!("B") => 0.1,
      o!("C") => 0.3,
      o!("D") => 0.1,
    };

    assert_eq!(expected.len(), actual.len());
    for (name, expected_div) in &expected {
      assert_abs_diff_eq!(actual[name], *expected_div, epsilon = 1e-8);
    }

    Ok(())
  }
}
