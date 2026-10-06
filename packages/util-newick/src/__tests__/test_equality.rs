#[cfg(test)]
mod tests {
  use crate::__tests__::test_parse_syntax::tests::helpers::{caterpillar, on_small_stack};
  use crate::parse::newick_from_string;
  use crate::types::{NewickEdgeData, NewickGraph, NewickHybrid, NewickNodeData, NewickReadOptions};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::reordered(            "((A,B),(C,D));",        "((D,C),(B,A));",        true)]
  #[case::other_topology(       "((A,B),C);",            "(A,(B,C));",            false)]
  #[case::other_length(         "(A:1,B:2);",            "(A:2,B:1);",            false)]
  #[case::other_annotation(     "(A[&x=1],B);",          "(A[&x=2],B);",          false)]
  #[case::same_multiset(        "((A,A),(A,B));",        "((A,B),(A,A));",        true)]
  #[case::different_multiset(   "((A,A),(B,B));",        "((A,B),(A,B));",        false)]
  #[trace]
  fn test_equality_unordered(#[case] left: &str, #[case] right: &str, #[case] expected: bool) {
    let options = NewickReadOptions::default();
    let left = newick_from_string(left, &options).unwrap();
    let right = newick_from_string(right, &options).unwrap();

    assert_eq!(expected, left == right);
  }

  #[test]
  fn test_equality_graph_with_cycle_is_reflexive() {
    let mut graph = NewickGraph::new();
    graph.root = graph.add_node(NewickNodeData::new());
    let hybrid = graph.add_node(NewickNodeData {
      hybrid: Some(NewickHybrid { kind: None, index: 1 }),
      ..NewickNodeData::new()
    });
    graph.add_edge(graph.root, hybrid, NewickEdgeData::new());
    graph.add_edge(hybrid, hybrid, NewickEdgeData::new());

    assert_eq!((true, false), (graph == graph.clone(), graph == NewickGraph::new()));
  }

  #[test]
  fn test_equality_deep_tree_on_small_stack() {
    let input = caterpillar(100_000);

    let equal = on_small_stack(move || {
      let graph = newick_from_string(&input, &NewickReadOptions::default()).unwrap();
      graph == graph.clone()
    });

    assert!(equal);
  }

  #[test]
  fn test_equality_stacked_hybrids_visit_each_node_once() {
    let input = (1..=40).fold("A".to_owned(), |inner, i| format!("(({inner})#H{i},(#H{i},B{i}))"));
    let options = NewickReadOptions { enewick: true };
    let graph = newick_from_string(&format!("{input};"), &options).unwrap();

    assert_eq!(graph, graph.clone());
  }
}
