#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::{read_error, read_with, summary};
  use crate::dialect::NewickDialect;
  use crate::model::comment::{EdgeComment, EdgeField, NewickComment, ValueSide};
  use crate::read::options::NewickReadOptions;
  use helpers::options;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::cardona_example(   "(A,B,((C,(Y)x#H1)c,(x#H1,D)d)e)f;",   vec!["f ", "A ", "B ", "e ", "c ", "C ", "x#H1  | ", "Y ", "d ", "D "])]
  #[case::acceptor(          "((A)x##LGT1,(x#LGT1,B));",            vec!["- ", "x#LGT1  | acceptor", "A ", "- ", "B "])]
  #[case::unnamed(           "((A)#H1,(#H1,B));",                   vec!["- ", "-#H1  | ", "A ", "- ", "B "])]
  #[case::later_name(        "((C)#H1,(X#H1,D));",                  vec!["- ", "X#H1  | ", "C ", "- ", "D "])]
  #[case::quoted_name(       "((C)'x y'#H1,('x y'#H1,D));",         vec!["- ", "x y#H1  | ", "C ", "- ", "D "])]
  #[case::quoted_is_name(    "('A#1',B);",                          vec!["- ", "A#1 ", "B "])]
  #[case::non_ascii_digits(  "(A#\u{663},B);",                      vec!["- ", "A#\u{663} ", "B "])]
  #[case::single_copy(       "(EPI_ISL#402124,B);",                 vec!["- ", "EPI_ISL#402124 ", "B "])]
  #[case::branch_lengths(    "((A)x#H1:1,(x#H1:2,B));",             vec!["- ", "x#H1 :2 | :1", "A ", "- ", "B "])]
  #[trace]
  fn test_read_networks_enewick(#[case] input: &str, #[case] expected: Vec<&str>) {
    assert_eq!(expected, summary(&read_with(input, &options(NewickDialect::ENEWICK)).graph));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::different_names(  "((C)x#H1,(y#H1,D));",     r#"line 1, column 11: The occurrences of the hybrid node #H1 have different names: "x" and "y""#)]
  #[case::children_twice(   "((A)#H1,(B)#H1);",        "line 1, column 12: The hybrid node #H1 has children in more than one of its occurrences")]
  #[case::parallel_edges(   "((C)#H1,#H1);",           "line 1, column 1: node 2 has more than one edge to node 1")]
  #[case::self_loop(        "(A,#H1)#H1;",             "line 1, column 1: The root, node 1, has a parent")]
  #[case::cycle(            "((#H2)#H1,(#H1)#H2);",    "line 1, column 1: The graph contains a cycle")]
  #[case::index_too_large(  "(A#H4294967296,B);",      r##"line 1, column 3: The hybrid node index in "#H4294967296" is not a valid index: number too large to fit in target type"##)]
  #[trace]
  fn test_read_networks_invalid(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(expected, read_error(input, &options(NewickDialect::ENEWICK)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::three_fields(       "(A:1:80:0.3,B::95)R:0.5;",            vec!["R :0.5", "A :1 support=80(field) p=0.3", "B support=95(field)"])]
  #[case::empty_fields(       "(A::,B:::0.5);",                      vec!["- ", "A ", "B p=0.5"])]
  #[case::label_and_field(    "((A,B)90:1:80,C);",                   vec!["- ", "90 :1 support=80(field)", "A ", "B ", "C "])]
  #[case::label_support(      "((A,B)90:1,C);",                      vec!["- ", "- :1 support=90(label)", "A ", "B ", "C "])]
  #[case::phylonet_gamma(     "((A,(B)#H1:::0.3)C,(#H1:::0.7,D)E)F;", vec!["F ", "C ", "A ", "-#H1 p=0.3 | p=0.7", "B ", "E ", "D "])]
  #[trace]
  fn test_read_networks_rich_fields(#[case] input: &str, #[case] expected: Vec<&str>) {
    assert_eq!(expected, summary(&read_with(input, &options(NewickDialect::RICH)).graph));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::weight_then_rooting(  "[&W 0.25][&U](A,B);",  (Some(false), Some(0.25)))]
  #[case::rooting_then_weight(  "[&r][&w 1](A,B);",     (Some(true),  Some(1.0)))]
  #[case::none(                 "(A,B);",               (None,        None))]
  #[trace]
  fn test_read_networks_rich_tree_comments(#[case] input: &str, #[case] expected: (Option<bool>, Option<f64>)) {
    let graph = read_with(input, &options(NewickDialect::RICH)).graph;

    assert_eq!(expected, (graph.rooted(), graph.weight()));
  }

  #[test]
  fn test_read_networks_rich_field_comments() {
    let tree = read_with("(A:1[c]:[d]80,B);", &options(NewickDialect::RICH));

    let edge = tree.graph.edge(tree.graph.child_edges(tree.graph.root())[0]).data();
    let expected = vec![
      EdgeComment::new(
        EdgeField::Length,
        ValueSide::AfterValue,
        NewickComment::Plain("c".to_owned()),
      ),
      EdgeComment::new(
        EdgeField::Support,
        ValueSide::BeforeValue,
        NewickComment::Plain("d".to_owned()),
      ),
    ];
    assert_eq!(expected, edge.comments().to_vec());
  }

  #[test]
  fn test_read_networks_second_rooting_comment_is_error() {
    assert_eq!(
      "line 1, column 5: The tree has more than one rooting comment",
      read_error("[&R][&U](A);", &options(NewickDialect::RICH))
    );
  }

  #[test]
  fn test_read_networks_hash_is_name_outside_network_dialects() {
    let actual = summary(&read_with("((A)x#H1,(x#H1,B));", &NewickReadOptions::default()).graph);

    assert_eq!(vec!["- ", "x#H1 ", "A ", "- ", "x#H1 ", "B "], actual);
  }

  mod helpers {
    use crate::dialect::NewickDialect;
    use crate::read::options::NewickReadOptions;

    pub(super) fn options(dialect: NewickDialect) -> NewickReadOptions {
      NewickReadOptions {
        dialect,
        ..NewickReadOptions::default()
      }
    }
  }
}
