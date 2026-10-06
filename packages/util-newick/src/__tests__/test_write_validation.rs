#[cfg(test)]
mod tests {
  use crate::__tests__::test_parse_syntax::tests::helpers::{caterpillar, on_small_stack};
  use crate::parse::newick_from_string;
  use crate::types::{
    NewickEdgeData, NewickGraph, NewickHybrid, NewickNodeData, NewickReadOptions, NewickValue, NewickWriteOptions,
    NwkStyle,
  };
  use crate::write::newick_to_string;
  use helpers::{options, star, write_error};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  #[rustfmt::skip]
  #[rstest]
  #[case::hash(           "A#1",    "('A#1',B);")]
  #[case::empty(          "",       "('',B);")]
  #[case::spaces(         "a b",    "('a b',B);")]
  #[case::apostrophe(     "it's",   "('it''s',B);")]
  #[case::numeric_leaf(   "123",    "(123,B);")]
  #[trace]
  fn test_write_validation_leaf_name_roundtrip(#[case] name: &str, #[case] expected: &str) {
    let g = star(NewickNodeData::new(), vec![NewickNodeData::new().with_name(name), NewickNodeData::new().with_name("B")]);

    let written = newick_to_string(&g, &options(NwkStyle::Plain)).unwrap();
    let parsed = newick_from_string(&written, &NewickReadOptions::default()).unwrap();

    assert_eq!((expected, true), (written.as_str(), g.eq_ordered(&parsed)));
  }

  #[test]
  fn test_write_validation_branch_annotation_without_length_stays_on_branch() {
    let mut edge = NewickEdgeData::new();
    edge.branch_attrs.insert("rate".to_owned(), NewickValue::Number(1.5));
    let mut g = star(NewickNodeData::new(), vec![NewickNodeData::new().with_name("A")]);
    g.edges[0].data = edge;

    let written = newick_to_string(&g, &options(NwkStyle::Beast)).unwrap();
    let parsed = newick_from_string(&written, &NewickReadOptions::default()).unwrap();

    assert_eq!(("(A:[&rate=1.5]);", true), (written.as_str(), g.eq_ordered(&parsed)));
  }

  #[test]
  fn test_write_validation_unnamed_leaf_annotation_stays_on_node() {
    let mut leaf = NewickNodeData::new();
    leaf.node_attrs.insert("a".to_owned(), NewickValue::Number(1.0));
    let mut g = star(NewickNodeData::new(), vec![leaf, NewickNodeData::new().with_name("B")]);
    g.edges[0].data.branch_length = Some(1.0);

    let written = newick_to_string(&g, &options(NwkStyle::Beast)).unwrap();
    let parsed = newick_from_string(&written, &NewickReadOptions::default()).unwrap();

    assert_eq!(("([&a=1]:1,B);", true), (written.as_str(), g.eq_ordered(&parsed)));
  }

  #[test]
  fn test_write_validation_hybrid_name_ending_in_hash_roundtrip() {
    let enewick = NewickReadOptions { enewick: true };
    let mut g = newick_from_string("(A,(B)x#H1,(x#H1,C));", &enewick).unwrap();
    let hybrid = g.nodes.iter().position(|node| node.hybrid.is_some()).unwrap();
    g.nodes[hybrid] = NewickNodeData {
      hybrid: Some(NewickHybrid {
        kind: Some("H".to_owned()),
        index: 1,
      }),
      children: g.nodes[hybrid].children.clone(),
      ..NewickNodeData::new().with_name("x#")
    };

    let written = newick_to_string(&g, &options(NwkStyle::Plain)).unwrap();
    let parsed = newick_from_string(&written, &enewick).unwrap();

    assert_eq!(
      ("(A,(B)'x#'#H1,('x#'#H1,C));", true),
      (written.as_str(), g.eq_ordered(&parsed))
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::leaf_support(        star(NewickNodeData::new(), vec![NewickNodeData::new().with_support(0.9)]),           "When writing Newick: When writing the label of node 1: A leaf cannot carry a support value, because Newick reads a leaf label as a name")]
  #[case::numeric_internal(    star(NewickNodeData::new().with_name("123"), vec![NewickNodeData::new()]),           "When writing Newick: When writing the label of node 0 ('123'): The internal node name \"123\" would be read back as a support value")]
  #[case::empty_graph(         NewickGraph::default(),                                                              "When writing Newick: The root is node 0, but the graph has 0 nodes")]
  #[trace]
  fn test_write_validation_rejects_unrepresentable_graph(#[case] graph: NewickGraph, #[case] expected: &str) {
    assert_eq!(expected, write_error(&graph, &options(NwkStyle::Plain)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::nan(       f64::NAN,       "When writing Newick: When writing the branch above node 1 ('A'): Newick cannot represent the number NaN")]
  #[case::infinite(  f64::INFINITY,  "When writing Newick: When writing the branch above node 1 ('A'): Newick cannot represent the number inf")]
  #[trace]
  fn test_write_validation_rejects_non_finite_branch_length(#[case] length: f64, #[case] expected: &str) {
    let mut g = star(NewickNodeData::new(), vec![NewickNodeData::new().with_name("A")]);
    g.edges[0].data.branch_length = Some(length);

    assert_eq!(expected, write_error(&g, &options(NwkStyle::Plain)));
  }

  #[test]
  fn test_write_validation_rejects_two_parents_without_hybrid_tag() {
    let mut g = star(
      NewickNodeData::new(),
      vec![NewickNodeData::new(), NewickNodeData::new().with_name("C")],
    );
    g.add_edge(1, 2, NewickEdgeData::new());

    assert_eq!(
      "When writing Newick: node 2 ('C') has 2 parents, but only a hybrid node can have more than one parent",
      write_error(&g, &options(NwkStyle::Plain))
    );
  }

  #[test]
  fn test_write_validation_rejects_zero_significant_digits() {
    let g = star(NewickNodeData::new(), vec![NewickNodeData::new().with_name("A")]);
    let zero_digits = NewickWriteOptions {
      significant_digits: Some(0),
      ..options(NwkStyle::Plain)
    };

    assert_eq!(
      "When writing Newick: The number of significant digits must be at least 1",
      write_error(&g, &zero_digits)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small_negative(  -0.00001,  "(A:-0);")]
  #[case::small_positive(  0.00049,   "(A:0);")]
  #[case::rounds_up(       0.0006,    "(A:0.001);")]
  #[trace]
  fn test_write_validation_decimal_digits_apply_to_small_numbers(#[case] length: f64, #[case] expected: &str) {
    let mut g = star(NewickNodeData::new(), vec![NewickNodeData::new().with_name("A")]);
    g.edges[0].data.branch_length = Some(length);
    let three_decimals = NewickWriteOptions {
      decimal_digits: Some(3),
      ..options(NwkStyle::Plain)
    };

    assert_eq!(expected, newick_to_string(&g, &three_decimals).unwrap());
  }

  #[test]
  fn test_write_validation_deep_tree_on_small_stack() {
    let input = caterpillar(100_000);
    let expected = format!("{};", input.strip_suffix(":0.01;").unwrap());

    let written = on_small_stack(move || {
      let g = newick_from_string(&input, &NewickReadOptions::default()).unwrap();
      newick_to_string(&g, &options(NwkStyle::Plain)).unwrap()
    });

    assert_eq!(expected, written);
  }

  #[test]
  fn test_write_validation_attributes_keep_written_order() {
    let mut leaf = NewickNodeData::new().with_name("A");
    leaf.node_attrs = BTreeMap::from([
      ("b".to_owned(), NewickValue::Boolean(true)),
      ("a".to_owned(), NewickValue::Boolean(false)),
    ]);
    let g = star(NewickNodeData::new(), vec![leaf]);

    assert_eq!(
      "(A[&a=FALSE,b=TRUE]);",
      newick_to_string(&g, &options(NwkStyle::Beast)).unwrap()
    );
  }

  mod helpers {
    use crate::types::{NewickEdgeData, NewickGraph, NewickNodeData, NewickWriteOptions, NwkStyle};
    use crate::write::newick_to_string;

    pub(super) fn options(style: NwkStyle) -> NewickWriteOptions {
      NewickWriteOptions {
        style,
        significant_digits: None,
        decimal_digits: None,
      }
    }

    pub(super) fn star(root: NewickNodeData, leaves: Vec<NewickNodeData>) -> NewickGraph {
      let mut graph = NewickGraph::new();
      graph.root = graph.add_node(root);
      for leaf in leaves {
        let child = graph.add_node(leaf);
        graph.add_edge(graph.root, child, NewickEdgeData::new());
      }
      graph
    }

    pub(super) fn write_error(graph: &NewickGraph, options: &NewickWriteOptions) -> String {
      format!("{:#}", newick_to_string(graph, options).unwrap_err())
    }
  }
}
