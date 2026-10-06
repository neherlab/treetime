#[cfg(test)]
pub(crate) mod tests {
  use crate::types::{NewickLabel, NewickValue};
  use helpers::{caterpillar, names, on_small_stack, parse, parse_error};
  use indoc::indoc;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;

  #[rustfmt::skip]
  #[rstest]
  #[case::missing_semicolon(  "(A:1,B:2)C",        "(A:1,B:2)C;")]
  #[case::leading_dot(        "(A:.5,B:2);",       "(A:0.5,B:2);")]
  #[case::plus_sign(          "(A:+1,B:2);",       "(A:1,B:2);")]
  #[case::byte_order_mark(    "\u{feff}(A,B);",    "(A,B);")]
  #[case::trailing_comment(   "(A,B);[end]",       "(A,B);")]
  #[case::spaces_around(      " ( A , B ) ; \n",   "(A,B);")]
  #[trace]
  fn test_parse_syntax_accepts_input_v0_accepts(#[case] input: &str, #[case] equivalent: &str) {
    let expected = parse(equivalent).unwrap();

    let actual = parse(input).unwrap();

    assert!(expected.eq_ordered(&actual), "{input:?} differs from {equivalent:?}");
  }

  #[test]
  fn test_parse_syntax_comment_before_open_belongs_to_that_node() {
    let g = parse("[c](A,B);").unwrap();

    assert_eq!(vec!["[c]"], g.nodes[g.root].raw_comments);
  }

  #[test]
  fn test_parse_syntax_comment_after_open_belongs_to_first_child() {
    let g = parse("([c]A,B);").unwrap();

    let a = g.nodes.iter().find(|node| node.name() == Some("A")).unwrap();
    assert_eq!(vec!["[c]"], a.raw_comments);
  }

  #[test]
  fn test_parse_syntax_comment_after_rooting_flag() {
    let g = parse("[&R][c](A,B);").unwrap();

    assert_eq!(
      (Some(true), vec!["[c]".to_owned()]),
      (g.rooted, g.nodes[g.root].raw_comments.clone())
    );
  }

  #[test]
  fn test_parse_syntax_lone_semicolon_is_one_unnamed_node() {
    let g = parse(";").unwrap();

    assert_eq!((1, None), (g.nodes.len(), g.nodes[g.root].label.clone()));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::close_bracket(  r#"(A[&a="x]y"],B);"#,  "x]y")]
  #[case::open_bracket(   r#"(A[&a="x[y"],B);"#,  "x[y")]
  #[case::semicolon(      r#"(A[&a="x;y"],B);"#,  "x;y")]
  #[trace]
  fn test_parse_syntax_brackets_inside_quoted_annotation(#[case] input: &str, #[case] expected: &str) {
    let g = parse(input).unwrap();

    let a = g.nodes.iter().find(|node| node.name() == Some("A")).unwrap();
    assert_eq!(Some(&NewickValue::String(expected.to_owned())), a.node_attrs.get("a"));
  }

  #[test]
  fn test_parse_syntax_unnamed_leaf_comment_before_colon_is_node_annotation() {
    let g = parse("([&a=1]:1,B);").unwrap();

    let edge = &g.edges[g.nodes[g.root].children[0]];
    let actual = (&g.nodes[edge.child].node_attrs, &edge.data.branch_attrs);
    assert_eq!(
      (
        &BTreeMap::from([("a".to_owned(), NewickValue::Number(1.0))]),
        &BTreeMap::new()
      ),
      actual
    );
  }

  #[test]
  fn test_parse_syntax_branch_annotation_without_length() {
    let g = parse("(A:[&rate=1.5],B);").unwrap();

    let edge = &g.edges[g.nodes[g.root].children[0]];
    let expected = (None, BTreeMap::from([("rate".to_owned(), NewickValue::Number(1.5))]));
    assert_eq!(expected, (edge.data.branch_length, edge.data.branch_attrs.clone()));
  }

  #[test]
  fn test_parse_syntax_hash_is_part_of_name_by_default() {
    let g = parse("(A#1,'B#2',EPI_ISL#402124);").unwrap();

    assert_eq!(vec!["A#1", "B#2", "EPI_ISL#402124"], names(&g));
  }

  #[test]
  fn test_parse_syntax_quoted_numeric_internal_label_is_support() {
    let g = parse("(A,B)'95';").unwrap();

    assert_eq!(Some(NewickLabel::Support(95.0)), g.nodes[g.root].label);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::extra_close(     "(A,B));",      "At line 1, column 6: unexpected ')' without a matching '('")]
  #[case::unclosed(        "((A,B);",      "At line 1, column 1: the '(' is never closed")]
  #[case::space_in_label(  "(A,B)C D;",    r#"At line 1, column 8: unexpected label "D": the node already has the label "C". Quote a label that contains spaces or punctuation"#)]
  #[case::label_before(    "A(B,C);",      "At line 1, column 2: unexpected '(': a node's '(' must come before its label and branch length")]
  #[case::two_lengths(     "(A:1:2,B);",   "At line 1, column 5: unexpected ':': the node already has a branch length")]
  #[case::label_after(     "(A:1B,C);",    r#"At line 1, column 5: unexpected label "B" after the branch length"#)]
  #[case::comma_outside(   "A,B;",         "At line 1, column 2: unexpected ',' outside of parentheses")]
  #[case::empty(           "",             "The input contains no tree")]
  #[trace]
  fn test_parse_syntax_error_names_position(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(format!("Failed to parse Newick string: {expected}"), parse_error(input));
  }

  #[test]
  fn test_parse_syntax_text_after_semicolon_is_error() {
    let expected = indoc! {"
      Failed to parse Newick string:  --> 1:7
        |
      1 | (A,B);(C,D);
        |       ^---
        |
        = expected end of input or comment"};

    assert_eq!(expected, parse_error("(A,B);(C,D);"));
  }

  #[test]
  fn test_parse_syntax_deep_tree_on_small_stack() {
    let input = caterpillar(100_000);

    let node_count = on_small_stack(move || parse(&input).unwrap().nodes.len());

    assert_eq!(200_001, node_count);
  }

  #[test]
  fn test_parse_syntax_deeply_nested_comment_on_small_stack() {
    let input = format!("(A{}{},B);", "[".repeat(100_000), "]".repeat(100_000));

    let raw_comment_length = on_small_stack(move || {
      let g = parse(&input).unwrap();
      g.nodes
        .iter()
        .map(|node| node.raw_comments.concat().len())
        .sum::<usize>()
    });

    assert_eq!(200_000, raw_comment_length);
  }

  pub(crate) mod helpers {
    use crate::parse::newick_from_string;
    use crate::types::{NewickGraph, NewickReadOptions};
    use std::fmt::Write;
    use std::thread;

    pub(super) fn parse(input: &str) -> eyre::Result<NewickGraph> {
      newick_from_string(input, &NewickReadOptions::default())
    }

    pub(super) fn parse_error(input: &str) -> String {
      format!("{:#}", parse(input).unwrap_err())
    }

    pub(super) fn names(graph: &NewickGraph) -> Vec<&str> {
      graph.nodes.iter().filter_map(|node| node.name()).collect()
    }

    pub(crate) fn caterpillar(depth: usize) -> String {
      let mut text = "(".repeat(depth);
      text.push_str("A:0.01");
      for i in 0..depth {
        write!(text, ",L{i}:0.01):0.01").unwrap();
      }
      text.push(';');
      text
    }

    pub(crate) fn on_small_stack<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
      thread::Builder::new()
        .stack_size(2 << 20)
        .spawn(f)
        .unwrap()
        .join()
        .unwrap()
    }
  }
}
