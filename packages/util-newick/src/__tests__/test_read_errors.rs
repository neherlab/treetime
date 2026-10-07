#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::{
    caterpillar, on_small_stack, read_error, read_with, summary,
  };
  use crate::dialect::NewickDialect;
  use crate::model::comment::NewickComment;
  use crate::read::error::NewickErrorKind;
  use crate::read::options::{NewickReadOptions, ReadMode};
  use crate::read::stream::newick_from_reader;
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::extra_close(       "(A,B));",       "line 1, column 6: unexpected ')' without a matching '('")]
  #[case::unclosed(          "((A,B);",       "line 1, column 1: the '(' is never closed")]
  #[case::space_in_label(    "(A,B)C D;",     "line 1, column 8: expected comment, ')', ',', ':', ';'")]
  #[case::label_before_open( "A(B,C);",       "line 1, column 2: expected comment, ')', ',', ':', ';'")]
  #[case::two_lengths(       "(A:1:2,B);",    "line 1, column 5: expected comment, ')', ',', ';'")]
  #[case::label_after_colon( "(A:1B,C);",     "line 1, column 5: expected comment, ')', ',', ';'")]
  #[case::comma_outside(     "A,B;",          "line 1, column 2: unexpected ',' outside of parentheses")]
  #[case::empty(             "",              "line 1, column 1: The input contains no tree")]
  #[case::only_comments(     " [c] ",         "line 1, column 1: The input contains no tree")]
  #[case::no_semicolon(      "(A,B)",         "line 1, column 6: The tree does not end with ';'")]
  #[case::second_tree(       "(A,B);(C,D);",  "line 1, column 7: The input contains more than one tree; read it tree by tree with newick_trees()")]
  #[case::text_after_tree(   "(A,B);xyz",     "line 1, column 7: The input contains more than one tree; read it tree by tree with newick_trees()")]
  #[case::second_line(       "(A,\nB C);",    "line 2, column 3: expected comment, ')', ',', ':', ';'")]
  #[case::byte_order_mark(   "\u{feff}(A);",  "line 1, column 1: unexpected byte order mark")]
  #[trace]
  fn test_read_errors_strict(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(expected, read_error(input, &NewickReadOptions::default()));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::missing_semicolon(  "(A:1,B:2)C",       vec!["C ", "A :1", "B :2"])]
  #[case::byte_order_mark(    "\u{feff}(A,B);",   vec!["- ", "A ", "B "])]
  #[case::trailing_comment(   "(A,B);[end]",      vec!["- ", "A ", "B "])]
  #[case::trailing_newlines(  "(A,B);\n\n",       vec!["- ", "A ", "B "])]
  #[trace]
  fn test_read_errors_tolerant_accepts(#[case] input: &str, #[case] expected: Vec<&str>) {
    let options = NewickReadOptions {
      mode: ReadMode::Tolerant,
      ..NewickReadOptions::default()
    };

    assert_eq!(expected, summary(&read_with(input, &options).graph));
  }

  #[test]
  fn test_read_errors_trailing_comment_strict() {
    assert_eq!(
      vec!["- ", "A ", "B "],
      summary(&read_with("(A,B);[end]", &NewickReadOptions::default()).graph)
    );
  }

  #[test]
  fn test_read_errors_invalid_utf8() {
    let actual = newick_from_reader(b"(A,\xff);".as_slice(), &NewickReadOptions::default()).unwrap_err();

    assert_eq!(
      (
        NewickErrorKind::InvalidUtf8,
        "line 1, column 4: The input is not valid UTF-8".to_owned()
      ),
      (actual.kind, actual.to_string())
    );
  }

  #[test]
  fn test_read_errors_syntax_error_of_the_given_dialect() {
    let options = NewickReadOptions {
      dialect: NewickDialect::NHX,
      ..NewickReadOptions::default()
    };

    let actual = read_error("(A[&a=1]B);", &options);

    assert_eq!(
      "line 1, column 9: expected comment, ')', ',', ':', ';', NHX annotation",
      actual
    );
  }

  #[test]
  fn test_read_errors_deep_tree_on_small_stack() {
    let input = caterpillar(100_000);

    let node_count = on_small_stack(move || read_with(&input, &NewickReadOptions::default()).graph.node_count());

    assert_eq!(200_001, node_count);
  }

  #[test]
  fn test_read_errors_deeply_nested_comment_on_small_stack() {
    let input = format!("(A{}{},B);", "[".repeat(100_000), "]".repeat(100_000));

    let length = on_small_stack(move || {
      let tree = read_with(&input, &NewickReadOptions::default());
      let leaf = tree.graph.children(tree.graph.root()).next().unwrap();
      match &tree.graph.node(leaf).comments()[0].comment {
        NewickComment::Plain(text) => text.len(),
        other => panic!("unexpected comment {other:?}"),
      }
    });

    assert_eq!(199_998, length);
  }

  #[test]
  fn test_read_errors_deeply_nested_array_on_small_stack() {
    let input = format!("(A[&a={}1{}],B);", "{".repeat(100_000), "}".repeat(100_000));
    let options = NewickReadOptions {
      dialect: NewickDialect::BEAST,
      ..NewickReadOptions::default()
    };

    let node_count = on_small_stack(move || read_with(&input, &options).graph.node_count());

    assert_eq!(3, node_count);
  }
}
