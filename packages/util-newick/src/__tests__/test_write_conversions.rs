#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::read_with;
  use crate::dialect::{NewickAnnotations, NewickDialect, NewickStructure};
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
  #[case::beast_comments(           DataKind::BeastComments,          "DKDD DKDD DKDD")]
  #[case::nhx_comments(             DataKind::NhxComments,            "DDKD DDKD DDKD")]
  #[case::mrbayes_comments(         DataKind::MrBayesComments,        "DDDK DDDK DDDK")]
  #[case::ampersand_comments(       DataKind::AmpersandPlainComments, "KFFF KFFF KFFF")]
  #[case::rooting_plain_comments(   DataKind::RootingPlainComments,   "KFKF KFKF FFFF")]
  #[case::weight_plain_comments(    DataKind::WeightPlainComments,    "KKKK KKKK FFFF")]
  #[case::rooting(                  DataKind::Rooting,                "DKDK DKDK KKKK")]
  #[case::weight(                   DataKind::Weight,                 "DDDD DDDD KKKK")]
  #[case::field_support(            DataKind::FieldSupport,           "DDDD DDDD KKKK")]
  #[case::probability(              DataKind::Probability,            "DDDD DDDD KKKK")]
  #[case::field_comments(           DataKind::FieldComments,          "DDDD DDDD KKKK")]
  #[case::hybrid_nodes(             DataKind::HybridNodes,            "FFFF KKKK KKKK")]
  #[trace]
  fn test_write_conversions_table(#[case] data: DataKind, #[case] expected: &str) {
    let actual: Vec<String> = NewickStructure::ALL
      .into_iter()
      .map(|structure| {
        NewickAnnotations::ALL
          .into_iter()
          .map(|annotations| match conversion(NewickDialect::new(structure, annotations), data) {
            Conversion::Keep => 'K',
            Conversion::Drop => 'D',
            Conversion::Fail => 'F',
          })
          .collect()
      })
      .collect();

    assert_eq!(expected, actual.join(" "));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beast_comments(     NewickDialect::BEAST,   "(A[&a=1][c]:1,B);",                 NewickDialect::CLASSIC, "(A[c]:1,B);")]
  #[case::nhx_comments(       NewickDialect::NHX,     "(A[&&NHX:S=x]:1,B);",               NewickDialect::BEAST,   "(A:1,B);")]
  #[case::mrbayes_comments(   NewickDialect::MRBAYES, "(A:1[&B r 0.1],B);",                NewickDialect::CLASSIC, "(A:1,B);")]
  #[case::rooting(            NewickDialect::BEAST,   "[&R](A,B);",                        NewickDialect::CLASSIC, "(A,B);")]
  #[case::weight(             NewickDialect::RICH,    "[&W 0.5](A,B);",                    NewickDialect::ENEWICK, "(A,B);")]
  #[case::field_support(      NewickDialect::RICH,    "(A:1:90,B);",                       NewickDialect::ENEWICK, "(A:1,B);")]
  #[case::probability(        NewickDialect::RICH,    "((A)#H1:1::0.3,(#H1:::0.7,B));",    NewickDialect::ENEWICK, "((A)#H1:1,(#H1,B));")]
  #[case::field_comments(     NewickDialect::RICH,    "(A:1:[c]90:[d],B);",                NewickDialect::CLASSIC, "(A:1,B);")]
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
  #[case::ampersand_comment(  NewickDialect::CLASSIC, "(A[&a=1],B);",      NewickDialect::BEAST, "When writing Newick: When writing node 0 ('A'): The comment Plain(\"&a=1\") cannot be written in the classic,beast dialect, which would read it as an annotation")]
  #[case::hybrid_nodes(       NewickDialect::ENEWICK, "((A)#H1,(#H1,B));", NewickDialect::NHX,   "When writing Newick: The classic,nhx dialect cannot hold hybrid nodes such as node 1; write networks with the enewick or rich structure")]
  #[case::rooting_comment(    NewickDialect::CLASSIC, "[&R]A;",            NewickDialect::RICH,  "When writing Newick: When writing node 0 ('A'): The comment [&R] at the start of the tree cannot be written in the rich,plain dialect, which would read it as the rooting comment")]
  #[case::unrooted_comment(   NewickDialect::CLASSIC, "[&U]A;",            NewickDialect::RICH,  "When writing Newick: When writing node 0 ('A'): The comment [&U] at the start of the tree cannot be written in the rich,plain dialect, which would read it as the rooting comment")]
  #[case::weight_comment(     NewickDialect::CLASSIC, "[&W 1]A;",          NewickDialect::RICH,  "When writing Newick: When writing node 0 ('A'): The comment [&W 1] at the start of the tree cannot be written in the rich,plain dialect, which would read it as the tree weight comment")]
  #[trace]
  fn test_write_conversions_fail(#[case] source: NewickDialect, #[case] input: &str, #[case] target: NewickDialect, #[case] expected: &str) {
    let actual = newick_to_string(&read_in(input, source), &NewickWriteOptions::new(target)).unwrap_err();

    assert_eq!(expected, format!("{actual:#}"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rooting_classic(      "[&R]A;",       NewickDialect::CLASSIC, "[&R]A;")]
  #[case::rooting_enewick(      "[&R]A;",       NewickDialect::ENEWICK, "[&R]A;")]
  #[case::weight_enewick(       "[&W 1]A;",     NewickDialect::ENEWICK, "[&W 1]A;")]
  #[case::after_children_rich(  "(A,B)[&R];",   NewickDialect::RICH,    "(A,B)[&R];")]
  #[trace]
  fn test_write_conversions_plain_tree_comment_kept(#[case] input: &str, #[case] target: NewickDialect, #[case] expected: &str) {
    let graph = read_in(input, NewickDialect::CLASSIC);

    let written = write_in(&graph, &NewickWriteOptions::new(target));

    assert_eq!((expected, true), (written.as_str(), graph.eq_ordered(&read_in(&written, target))));
  }

  #[test]
  fn test_write_conversions_network_occurrence_annotations() {
    let input = "((A[&segments={0,1}]:1,(B[&segments={0,1}]:1)#H0[&segments={0}]:0.5)[&segments={0,1}]:1,(#H0[&segments={1}]:0.7,C[&segments={0,1}]:1)[&segments={0,1}]:1);";
    let graph = read_in(input, NewickDialect::ENEWICK_BEAST);

    let written = write_in(&graph, &NewickWriteOptions::new(NewickDialect::ENEWICK_BEAST));

    assert_eq!(input, written);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::root_edge(    NewickWriteOptions { root_edge: false, ..NewickWriteOptions::default() },                                                  "(A:0.123456,B)95:0.5;",  "(A:0.123456,B);")]
  #[case::digits(       NewickWriteOptions { numbers: NumberFormat { significant_digits: Some(2), ..NumberFormat::default() }, ..NewickWriteOptions::default() }, "(A:0.123456,B)95:0.5;", "(A:0.12,B)95:0.5;")]
  #[case::underscores(  NewickWriteOptions { spaces: Spaces::Underscore, ..NewickWriteOptions::default() },                                         "('a b':1,c_d);",         "(a_b:1,c_d);")]
  #[trace]
  fn test_write_conversions_options(#[case] options: NewickWriteOptions, #[case] input: &str, #[case] expected: &str) {
    let written = write_in(&read_in(input, NewickDialect::CLASSIC), &options);

    assert_eq!(expected, written);
  }

  #[test]
  fn test_write_conversions_underscores_read_back_as_spaces() {
    let options = NewickWriteOptions {
      spaces: Spaces::Underscore,
      ..NewickWriteOptions::default()
    };
    let graph = read_in("('a b':1,'c d');", NewickDialect::CLASSIC);

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
        dialect,
        ..NewickReadOptions::default()
      };
      read_with(text, &options).graph
    }

    pub(super) fn write_in(graph: &NewickGraph, options: &NewickWriteOptions) -> String {
      newick_to_string(graph, options).unwrap()
    }
  }
}
