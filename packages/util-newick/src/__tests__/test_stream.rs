#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::summary;
  use crate::dialect::NewickDialect;
  use crate::model::comment::NewickComment;
  use crate::read::options::{NewickReadOptions, ReadMode};
  use crate::read::stream::newick_trees;
  use helpers::{names_and_errors, split_at_first_read};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[test]
  fn test_stream_reads_tree_by_tree() {
    let input = "(A,B);\n[between] (C,D);\n\n(E);\n[trailing]\n";

    let actual = names_and_errors(input.as_bytes(), &NewickReadOptions::default());

    assert_eq!(
      vec![
        Ok(vec!["- ".to_owned(), "A ".to_owned(), "B ".to_owned()]),
        Ok(vec!["- ".to_owned(), "C ".to_owned(), "D ".to_owned()]),
        Ok(vec!["- ".to_owned(), "E ".to_owned()]),
      ],
      actual
    );
  }

  #[test]
  fn test_stream_comment_between_trees_belongs_to_next_root() {
    let trees: Vec<_> = newick_trees(b"(A);[c](B);".as_slice(), NewickReadOptions::default())
      .map(Result::unwrap)
      .collect();

    let root = trees[1].graph.root();
    assert_eq!(
      vec![NewickComment::Plain("c".to_owned())],
      trees[1]
        .graph
        .node(root)
        .comments()
        .iter()
        .map(|comment| comment.comment.clone())
        .collect::<Vec<_>>()
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::strict(    ReadMode::Strict,   vec![Ok(vec!["- ".to_owned(), "A ".to_owned()]), Err("line 2, column 4: The tree does not end with ';'".to_owned())])]
  #[case::tolerant(  ReadMode::Tolerant, vec![Ok(vec!["- ".to_owned(), "A ".to_owned()]), Ok(vec!["- ".to_owned(), "B ".to_owned()])])]
  #[trace]
  fn test_stream_last_tree_without_semicolon(#[case] mode: ReadMode, #[case] expected: Vec<Result<Vec<String>, String>>) {
    let options = NewickReadOptions {
      mode,
      ..NewickReadOptions::default()
    };

    assert_eq!(expected, names_and_errors(b"(A);\n(B)\n".as_slice(), &options));
  }

  #[test]
  fn test_stream_stops_at_invalid_utf8() {
    let actual = names_and_errors(b"(A);\n(B\xff);(C);".as_slice(), &NewickReadOptions::default());

    assert_eq!(
      vec![
        Ok(vec!["- ".to_owned(), "A ".to_owned()]),
        Err("line 2, column 3: The input is not valid UTF-8".to_owned()),
      ],
      actual
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::comment_with_semicolon(  NewickDialect::Classic, "(C[x;y],D);",          5, vec!["- ", "C ", "D "])]
  #[case::quoted_semicolon(        NewickDialect::Classic, "('q;r',D);",           4, vec!["- ", "q;r ", "D "])]
  #[case::string_with_bracket(     NewickDialect::Beast,   "(C[&a=\"x];y\"],D);",  10, vec!["- ", "C ", "D "])]
  #[case::two_byte_character(      NewickDialect::Classic, "(C\u{e9},D);",         3, vec!["- ", "C\u{e9} ", "D "])]
  #[case::four_byte_character(     NewickDialect::Classic, "(C\u{1f600},D);",      4, vec!["- ", "C\u{1f600} ", "D "])]
  #[trace]
  fn test_stream_construct_across_read_boundary(
    #[case] dialect: NewickDialect,
    #[case] second_tree: &str,
    #[case] split: usize,
    #[case] expected: Vec<&str>,
  ) {
    let options = NewickReadOptions {
      dialects: vec![dialect],
      ..NewickReadOptions::default()
    };
    let input = split_at_first_read(second_tree, split);

    let actual = names_and_errors(input.as_bytes(), &options);

    let expected: Vec<String> = expected.into_iter().map(str::to_owned).collect();
    assert_eq!((2, Some(Ok(expected))), (actual.len(), actual.get(1).cloned()));
  }

  mod helpers {
    use super::summary;
    use crate::read::options::NewickReadOptions;
    use crate::read::stream::newick_trees;
    use std::io::Read;

    const FIRST_READ_SIZE: usize = 64 * 1024;

    pub(super) fn names_and_errors(reader: impl Read, options: &NewickReadOptions) -> Vec<Result<Vec<String>, String>> {
      newick_trees(reader, options.clone())
        .map(|tree| tree.map(|tree| summary(&tree.graph)).map_err(|error| error.to_string()))
        .collect()
    }

    pub(super) fn split_at_first_read(second_tree: &str, split: usize) -> String {
      let filler = "A".repeat(FIRST_READ_SIZE - split - "(,B);".len());
      format!("({filler},B);{second_tree}")
    }
  }
}
