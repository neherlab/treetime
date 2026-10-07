#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::{read_error, read_with, summary};
  use crate::dialect::NewickDialect;
  use crate::model::comment::{
    EdgeComment, EdgeField, LabelSide, MrBayesComment, MrBayesKind, NewickComment, NodeComment, ValueSide,
  };
  use crate::model::value::NewickValue;
  use crate::read::error::NewickWarning;
  use crate::read::options::{NewickReadOptions, ReadMode};
  use helpers::{first_leaf_comments, options, pairs};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::shortest_number(   "[&a=0.95]",                 NewickValue::Number(0.95))]
  #[case::integer(           "[&a=12]",                   NewickValue::Number(12.0))]
  #[case::negative_zero(     "[&a=-0]",                   NewickValue::Number(-0.0))]
  #[case::leading_zero(      "[&a=0123]",                 NewickValue::NumberText("0123".to_owned()))]
  #[case::trailing_zero(     "[&a=2020.50]",              NewickValue::NumberText("2020.50".to_owned()))]
  #[case::exponent(          "[&a=1e5]",                  NewickValue::NumberText("1e5".to_owned()))]
  #[case::big_integer(       "[&a=12345678901234567890]", NewickValue::NumberText("12345678901234567890".to_owned()))]
  #[case::overflow(          "[&a=1e400]",                NewickValue::NumberText("1e400".to_owned()))]
  #[case::infinity_word(     "[&a=-infinity]",            NewickValue::String("-infinity".to_owned()))]
  #[case::single_quoted(     "[&a='x,y']",                NewickValue::String("x,y".to_owned()))]
  #[case::double_quoted(     "[&a=\"x,y\"]",              NewickValue::String("x,y".to_owned()))]
  #[case::doubled_quote(     "[&a=\"say \"\"hi\"\"\"]",   NewickValue::String("say \"hi\"".to_owned()))]
  #[case::quoted_true(       "[&a=\"TRUE\"]",             NewickValue::String("TRUE".to_owned()))]
  #[case::bracket_in_string( "[&a=\"x]y\"]",              NewickValue::String("x]y".to_owned()))]
  #[case::apostrophe_word(   "[&a=Cote d'Ivoire]",        NewickValue::String("Cote d'Ivoire".to_owned()))]
  #[case::boolean(           "[&a=True]",                 NewickValue::Boolean(true))]
  #[case::boolean_false(     "[&a=FALSE]",                NewickValue::Boolean(false))]
  #[case::bare_key(          "[&a]",                      NewickValue::Boolean(true))]
  #[case::color(             "[&a=#Ff0080]",              NewickValue::Color([255, 0, 128]))]
  #[case::empty_array(       "[&a={}]",                   NewickValue::Array(vec![].into()))]
  #[case::quoted_brace(      "[&a={\"x}\",y}]",           NewickValue::Array(vec![NewickValue::String("x}".to_owned()), NewickValue::String("y".to_owned())].into()))]
  #[case::nested_arrays(     "[&a={{1,2},{}}]",           NewickValue::Array(vec![NewickValue::Array(vec![NewickValue::Number(1.0), NewickValue::Number(2.0)].into()), NewickValue::Array(vec![].into())].into()))]
  #[trace]
  fn test_read_comments_beast_value(#[case] comment: &str, #[case] expected: NewickValue) {
    let input = format!("(A{comment},B);");

    let actual = first_leaf_comments(&input, &options(NewickDialect::BEAST));

    assert_eq!(vec![NodeComment::new(LabelSide::AfterLabel, NewickComment::Beast(pairs(vec![("a", expected)])))], actual);
  }

  #[test]
  fn test_read_comments_beast_pairs_keep_order_and_repeated_keys() {
    let input = "(A[&b=1,\"a b\"=2,b=3],B);";

    let actual = first_leaf_comments(input, &options(NewickDialect::BEAST));

    let expected = NewickComment::Beast(pairs(vec![
      ("b", NewickValue::Number(1.0)),
      ("a b", NewickValue::Number(2.0)),
      ("b", NewickValue::Number(3.0)),
    ]));
    assert_eq!(vec![NodeComment::new(LabelSide::AfterLabel, expected)], actual);
  }

  #[test]
  fn test_read_comments_annotation_accessor_returns_last_value() {
    let tree = read_with("(A[&b=1][&b=3],B);", &options(NewickDialect::BEAST));
    let leaf = tree.graph.children(tree.graph.root()).next().unwrap();

    assert_eq!(Some(&NewickValue::Number(3.0)), tree.graph.node(leaf).annotation("b"));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::standard_tags(  "[&&NHX:S=human:T=9606:B=90:D=Y]",  vec![("S", NewickValue::String("human".to_owned())), ("T", NewickValue::Number(9606.0)), ("B", NewickValue::Number(90.0)), ("D", NewickValue::String("Y".to_owned()))])]
  #[case::color(          "[&&NHX:C=255.0.12]",               vec![("C", NewickValue::Color([255, 0, 12]))])]
  #[case::compound(       "[&&NHX:Ev=1>1>0>dup>x]",           vec![("Ev", NewickValue::Array(vec![NewickValue::String("1".to_owned()), NewickValue::String("1".to_owned()), NewickValue::String("0".to_owned()), NewickValue::String("dup".to_owned()), NewickValue::String("x".to_owned())].into()))])]
  #[case::custom_number(  "[&&NHX:date=2020.50]",             vec![("date", NewickValue::String("2020.50".to_owned()))])]
  #[case::bare_tag(       "[&&nhx:D:S=x]",                    vec![("D", NewickValue::Boolean(true)), ("S", NewickValue::String("x".to_owned()))])]
  #[case::comma_in_value( "[&&NHX:m=A1T,C2G]",                vec![("m", NewickValue::String("A1T,C2G".to_owned()))])]
  #[case::empty_value(    "[&&NHX:S=]",                       vec![("S", NewickValue::String(String::new()))])]
  #[case::trailing_colon( "[&&NHX:S=x:]",                     vec![("S", NewickValue::String("x".to_owned()))])]
  #[trace]
  fn test_read_comments_nhx_tags(#[case] comment: &str, #[case] expected: Vec<(&str, NewickValue)>) {
    let input = format!("(A{comment},B);");

    let actual = first_leaf_comments(&input, &options(NewickDialect::NHX));

    assert_eq!(vec![NodeComment::new(LabelSide::AfterLabel, NewickComment::Nhx(pairs(expected)))], actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::decimal(      "[&&NHX:B=high]",    r#"line 1, column 10: The NHX tag B needs a decimal number, but its value is "high""#)]
  #[case::integer(      "[&&NHX:T=1.5]",     r#"line 1, column 10: The NHX tag T needs an integer, but its value is "1.5""#)]
  #[case::duplication(  "[&&NHX:D=maybe]",   r#"line 1, column 10: The NHX tag D needs one of T, F, Y, N or ?, but its value is "maybe""#)]
  #[case::color_range(  "[&&NHX:C=300.0.0]", r#"line 1, column 10: The NHX tag C needs a color written as red.green.blue, but its value is "300.0.0""#)]
  #[case::parts(        "[&&NHX:B=1>2]",     r#"line 1, column 10: The NHX tag B needs a decimal number, but its value is "1>2""#)]
  #[trace]
  fn test_read_comments_nhx_type_mismatch_strict(#[case] comment: &str, #[case] expected: &str) {
    let input = format!("(A{comment},B);");

    assert_eq!(expected, read_error(&input, &options(NewickDialect::NHX)));
  }

  #[test]
  fn test_read_comments_nhx_type_mismatch_tolerant() {
    let options = NewickReadOptions {
      mode: ReadMode::Tolerant,
      ..options(NewickDialect::NHX)
    };

    let tree = read_with("(A[&&NHX:B=high],B);", &options);

    let leaf = tree.graph.children(tree.graph.root()).next().unwrap();
    let expected_comment = NewickComment::Nhx(pairs(vec![("B", NewickValue::String("high".to_owned()))]));
    let expected_warning = NewickWarning {
      offset: 9,
      line: 1,
      column: 10,
      message: r#"The NHX tag B needs a decimal number, but its value is "high""#.to_owned(),
    };
    assert_eq!(
      (
        vec![NodeComment::new(LabelSide::AfterLabel, expected_comment)],
        vec![expected_warning]
      ),
      (tree.graph.node(leaf).comments().to_vec(), tree.warnings)
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::branch(  "[&B TK02Brlens 0.1]",   MrBayesKind::B, "TK02Brlens", vec!["0.1"])]
  #[case::event(   "[&E ibr 2: 0.1 0.2]",   MrBayesKind::E, "ibr",        vec!["2:", "0.1", "0.2"])]
  #[case::node(    "[&N height]",           MrBayesKind::N, "height",     vec![])]
  #[trace]
  fn test_read_comments_mrbayes(#[case] comment: &str, #[case] kind: MrBayesKind, #[case] name: &str, #[case] values: Vec<&str>) {
    let input = format!("(A{comment},B);");

    let actual = first_leaf_comments(&input, &options(NewickDialect::MRBAYES));

    let expected = NewickComment::MrBayesMcmc(MrBayesComment {
      kind,
      name: name.to_owned(),
      values: values.into_iter().map(str::to_owned).collect(),
    });
    assert_eq!(vec![NodeComment::new(LabelSide::AfterLabel, expected)], actual);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beast_missing_value(   NewickDialect::BEAST,   "(A[&a=],B);",       "line 1, column 3: The annotation [&a=] does not follow the annotation syntax of the dialect")]
  #[case::beast_unclosed_string( NewickDialect::BEAST,   "(A[&a=\"x],B);",    "line 1, column 3: The annotation [&a=\"x] does not follow the annotation syntax of the dialect")]
  #[case::nhx_with_beast(        NewickDialect::NHX,     "(A[&a=1],B);",      "line 1, column 3: The annotation [&a=1] does not follow the annotation syntax of the dialect")]
  #[case::enewick_beast_malformed(NewickDialect::ENEWICK_BEAST, "(A[&a=],B);", "line 1, column 3: The annotation [&a=] does not follow the annotation syntax of the dialect")]
  #[case::mrbayes_with_beast(    NewickDialect::MRBAYES, "(A[&prob=1],B);",   "line 1, column 3: The annotation [&prob=1] does not follow the annotation syntax of the dialect")]
  #[trace]
  fn test_read_comments_malformed_annotation_strict(#[case] dialect: NewickDialect, #[case] input: &str, #[case] expected: &str) {
    assert_eq!(expected, read_error(input, &options(dialect)));
  }

  #[test]
  fn test_read_comments_malformed_annotation_tolerant() {
    let options = NewickReadOptions {
      mode: ReadMode::Tolerant,
      ..options(NewickDialect::BEAST)
    };

    let tree = read_with("(A[&a=],B);", &options);

    let leaf = tree.graph.children(tree.graph.root()).next().unwrap();
    assert_eq!(
      (
        vec![NodeComment::new(
          LabelSide::AfterLabel,
          NewickComment::Plain("&a=".to_owned())
        )],
        vec!["line 1, column 3: The annotation [&a=] does not follow the annotation syntax of the dialect".to_owned()]
      ),
      (
        tree.graph.node(leaf).comments().to_vec(),
        tree.warnings.iter().map(ToString::to_string).collect::<Vec<_>>()
      )
    );
  }

  #[test]
  fn test_read_comments_quotes_in_plain_comments_are_text() {
    let tree = read_with("(A[5\" tall],B[say \"hi\"]);", &NewickReadOptions::default());

    let comments: Vec<Vec<NodeComment>> = tree
      .graph
      .children(tree.graph.root())
      .map(|child| tree.graph.node(child).comments().to_vec())
      .collect();
    assert_eq!(
      vec![
        vec![NodeComment::new(
          LabelSide::AfterLabel,
          NewickComment::Plain("5\" tall".to_owned())
        )],
        vec![NodeComment::new(
          LabelSide::AfterLabel,
          NewickComment::Plain("say \"hi\"".to_owned())
        )],
      ],
      comments
    );
  }

  #[test]
  fn test_read_comments_positions() {
    let tree = read_with("[r]([a]A[b]:[c]1[d],[e]:2)Y[f];", &options(NewickDialect::CLASSIC));
    let graph = &tree.graph;
    let children: Vec<usize> = graph.children(graph.root()).collect();
    let edges: Vec<usize> = graph.child_edges(graph.root()).to_vec();

    let plain = |text: &str| NewickComment::Plain(text.to_owned());
    let actual = (
      graph.node(graph.root()).comments().to_vec(),
      graph.node(children[0]).comments().to_vec(),
      graph.edge(edges[0]).data().comments().to_vec(),
      graph.node(children[1]).comments().to_vec(),
    );
    let expected = (
      vec![
        NodeComment::new(LabelSide::BeforeLabel, plain("r")),
        NodeComment::new(LabelSide::AfterLabel, plain("f")),
      ],
      vec![
        NodeComment::new(LabelSide::BeforeLabel, plain("a")),
        NodeComment::new(LabelSide::AfterLabel, plain("b")),
      ],
      vec![
        EdgeComment::new(EdgeField::Length, ValueSide::BeforeValue, plain("c")),
        EdgeComment::new(EdgeField::Length, ValueSide::AfterValue, plain("d")),
      ],
      vec![NodeComment::new(LabelSide::AfterLabel, plain("e"))],
    );
    assert_eq!(expected, actual);
  }

  #[test]
  fn test_read_comments_branch_annotation_without_length() {
    let tree = read_with("(A:[&rate=1.5],B);", &options(NewickDialect::BEAST));

    let edge = tree.graph.edge(tree.graph.child_edges(tree.graph.root())[0]).data();
    let expected = vec![EdgeComment::new(
      EdgeField::Length,
      ValueSide::BeforeValue,
      NewickComment::Beast(pairs(vec![("rate", NewickValue::Number(1.5))])),
    )];
    assert_eq!((None, expected), (edge.branch_length(), edge.comments().to_vec()));
  }

  #[test]
  fn test_read_comments_rooting_comment_only_in_rooting_dialects() {
    let classic = read_with("[&R](A,B);", &options(NewickDialect::CLASSIC));
    let beast = read_with("[&R](A,B);", &options(NewickDialect::BEAST));

    let expected_classic = vec![NodeComment::new(
      LabelSide::AfterLabel,
      NewickComment::Plain("&R".to_owned()),
    )];
    assert_eq!(
      ((None, expected_classic), (Some(true), vec![])),
      (
        (
          classic.graph.rooted(),
          classic.graph.node(classic.graph.root()).comments().to_vec()
        ),
        (
          beast.graph.rooted(),
          beast.graph.node(beast.graph.root()).comments().to_vec()
        )
      )
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::enewick_plain(  NewickDialect::ENEWICK,       NewickComment::Plain("&a=1".to_owned()))]
  #[case::enewick_beast(  NewickDialect::ENEWICK_BEAST, NewickComment::Beast(pairs(vec![("a", NewickValue::Number(1.0))])))]
  #[trace]
  fn test_read_comments_enewick_annotation_convention(#[case] dialect: NewickDialect, #[case] expected: NewickComment) {
    let actual = first_leaf_comments("(A[&a=1],B);", &options(dialect));

    assert_eq!(vec![NodeComment::new(LabelSide::AfterLabel, expected)], actual);
  }

  #[test]
  fn test_read_comments_classic_keeps_beast_annotations_as_text() {
    let tree = read_with("(A[&a=1]:1,B);", &NewickReadOptions::default());

    assert_eq!(vec!["- ", "A :1", "B "], summary(&tree.graph));
  }

  mod helpers {
    use crate::__tests__::test_read_basic::tests::helpers::read_with;
    use crate::dialect::NewickDialect;
    use crate::model::comment::NodeComment;
    use crate::model::value::NewickValue;
    use crate::read::options::NewickReadOptions;

    pub(super) fn options(dialect: NewickDialect) -> NewickReadOptions {
      NewickReadOptions {
        dialect,
        ..NewickReadOptions::default()
      }
    }

    pub(super) fn pairs(pairs: Vec<(&str, NewickValue)>) -> Vec<(String, NewickValue)> {
      pairs.into_iter().map(|(key, value)| (key.to_owned(), value)).collect()
    }

    pub(super) fn first_leaf_comments(input: &str, options: &NewickReadOptions) -> Vec<NodeComment> {
      let tree = read_with(input, options);
      let leaf = tree.graph.children(tree.graph.root()).next().unwrap();
      tree.graph.node(leaf).comments().to_vec()
    }
  }
}
