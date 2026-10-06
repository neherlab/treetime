#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::read_with;
  use crate::dialect::NewickDialect;
  use crate::number::NumberFormat;
  use crate::read::options::NewickReadOptions;
  use crate::write::conversions::{Conversion, DataKind, conversion};
  use crate::write::newick::newick_to_string;
  use crate::write::options::{NewickWriteOptions, Spaces};
  use helpers::{read_in, write_in};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::beast_comments(     DataKind::BeastComments,          "DKDDDD")]
  #[case::nhx_comments(       DataKind::NhxComments,            "DDDKDD")]
  #[case::mrbayes_comments(   DataKind::MrBayesComments,        "DDKDDD")]
  #[case::ampersand_comments( DataKind::AmpersandPlainComments, "KFFFFF")]
  #[case::rooting(            DataKind::Rooting,                "DKKDDK")]
  #[case::weight(             DataKind::Weight,                 "DDDDDK")]
  #[case::field_support(      DataKind::FieldSupport,           "DDDDDK")]
  #[case::probability(        DataKind::Probability,            "DDDDDK")]
  #[case::field_comments(     DataKind::FieldComments,          "DDDDDK")]
  #[case::hybrid_nodes(       DataKind::HybridNodes,            "FFFFKK")]
  #[trace]
  fn test_write_conversions_table(#[case] data: DataKind, #[case] expected: &str) {
    let dialects = [
      NewickDialect::Classic,
      NewickDialect::Beast,
      NewickDialect::MrBayes,
      NewickDialect::Nhx,
      NewickDialect::ENewick,
      NewickDialect::Rich,
    ];

    let actual: String = dialects
      .iter()
      .map(|&dialect| match conversion(dialect, data) {
        Conversion::Keep => 'K',
        Conversion::Drop => 'D',
        Conversion::Fail => 'F',
      })
      .collect();

    assert_eq!(expected, actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beast_comments(     NewickDialect::Beast,   "(A[&a=1][c]:1,B);",                 NewickDialect::Classic, "(A[c]:1,B);")]
  #[case::nhx_comments(       NewickDialect::Nhx,     "(A[&&NHX:S=x]:1,B);",               NewickDialect::Beast,   "(A:1,B);")]
  #[case::mrbayes_comments(   NewickDialect::MrBayes, "(A:1[&B r 0.1],B);",                NewickDialect::Classic, "(A:1,B);")]
  #[case::rooting(            NewickDialect::Beast,   "[&R](A,B);",                        NewickDialect::Classic, "(A,B);")]
  #[case::weight(             NewickDialect::Rich,    "[&W 0.5](A,B);",                    NewickDialect::ENewick, "(A,B);")]
  #[case::field_support(      NewickDialect::Rich,    "(A:1:90,B);",                       NewickDialect::ENewick, "(A:1,B);")]
  #[case::probability(        NewickDialect::Rich,    "((A)#H1:1::0.3,(#H1:::0.7,B));",    NewickDialect::ENewick, "((A)#H1:1,(#H1,B));")]
  #[case::field_comments(     NewickDialect::Rich,    "(A:1:[c]90:[d],B);",                NewickDialect::Classic, "(A:1,B);")]
  #[trace]
  fn test_write_conversions_drop_only_that_data(
    #[case] source: NewickDialect,
    #[case] input: &str,
    #[case] target: NewickDialect,
    #[case] expected: &str,
  ) {
    let written = write_in(&read_in(input, source), &NewickWriteOptions::new(target));

    assert_eq!((expected, true), (written.as_str(), read_in(expected, source).eq_ordered(&read_in(&written, source))));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::ampersand_comment(  NewickDialect::Classic, "(A[&a=1],B);",          NewickDialect::Beast, "When writing Newick: When writing node 0 ('A'): The comment Plain(\"&a=1\") cannot be written in the beast dialect, which would read it as an annotation")]
  #[case::hybrid_nodes(       NewickDialect::ENewick, "((A)#H1,(#H1,B));",     NewickDialect::Nhx,   "When writing Newick: The nhx dialect cannot hold hybrid nodes such as node 1; write networks in the enewick or rich dialect")]
  #[trace]
  fn test_write_conversions_fail(#[case] source: NewickDialect, #[case] input: &str, #[case] target: NewickDialect, #[case] expected: &str) {
    let actual = newick_to_string(&read_in(input, source), &NewickWriteOptions::new(target)).unwrap_err();

    assert_eq!(expected, format!("{actual:#}"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::root_edge(    NewickWriteOptions { root_edge: false, ..NewickWriteOptions::default() },                                                  "(A:0.123456,B)95:0.5;",  "(A:0.123456,B);")]
  #[case::digits(       NewickWriteOptions { numbers: NumberFormat { significant_digits: Some(2), ..NumberFormat::default() }, ..NewickWriteOptions::default() }, "(A:0.123456,B)95:0.5;", "(A:0.12,B)95:0.5;")]
  #[case::underscores(  NewickWriteOptions { spaces: Spaces::Underscore, ..NewickWriteOptions::default() },                                         "('a b':1,c_d);",         "(a_b:1,c_d);")]
  #[trace]
  fn test_write_conversions_options(#[case] options: NewickWriteOptions, #[case] input: &str, #[case] expected: &str) {
    let written = write_in(&read_in(input, NewickDialect::Classic), &options);

    assert_eq!(expected, written);
  }

  #[test]
  fn test_write_conversions_underscores_read_back_as_spaces() {
    let options = NewickWriteOptions {
      spaces: Spaces::Underscore,
      ..NewickWriteOptions::default()
    };
    let graph = read_in("('a b':1,'c d');", NewickDialect::Classic);

    let written = newick_to_string(&graph, &options).unwrap();
    let reread = read_with(
      &written,
      &NewickReadOptions {
        underscores_as_spaces: true,
        ..NewickReadOptions::default()
      },
    );

    assert!(graph.eq_ordered(&reread.graph), "{written} reads back differently");
  }

  mod helpers {
    use crate::__tests__::test_read_basic::tests::helpers::read_with;
    use crate::dialect::NewickDialect;
    use crate::model::graph::NewickGraph;
    use crate::read::options::NewickReadOptions;
    use crate::write::newick::newick_to_string;
    use crate::write::options::NewickWriteOptions;

    pub(super) fn read_in(text: &str, dialect: NewickDialect) -> NewickGraph {
      let options = NewickReadOptions {
        dialects: vec![dialect],
        ..NewickReadOptions::default()
      };
      read_with(text, &options).graph
    }

    pub(super) fn write_in(graph: &NewickGraph, options: &NewickWriteOptions) -> String {
      newick_to_string(graph, options).unwrap()
    }
  }
}
