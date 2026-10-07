#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::{caterpillar, on_small_stack, read_with};
  use crate::dialect::{NewickAnnotations, NewickDialect, NewickStructure};
  use crate::model::comment::NewickComment;
  use crate::model::data::{NewickEdgeData, NewickNodeData, SupportSource};
  use crate::model::graph::NewickGraph;
  use crate::model::value::NewickValue;
  use crate::number::NumberFormat;
  use crate::read::options::NewickReadOptions;
  use crate::write::newick::{newick_to_string, write_newick_trees};
  use crate::write::options::{BranchAnnotations, NewickWriteOptions, Quoting, Spaces, SupportPlacement};
  use helpers::{
    commented_leaf, named_internal_with_support, named_root_support, reread, rewrite, star, star_with_edge, write_error,
  };
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::classic(           NewickDialect::CLASSIC, "(A:0.1,B:0.2)root;")]
  #[case::empty_names(       NewickDialect::CLASSIC, "(,);")]
  #[case::root_length(       NewickDialect::CLASSIC, "(A,B):0.5;")]
  #[case::support(           NewickDialect::CLASSIC, "((A,B)95:1,C)0.5;")]
  #[case::multi_support(     NewickDialect::CLASSIC, "((A,B)80.5/95,C);")]
  #[case::exponent(          NewickDialect::CLASSIC, "(A:1.0e-10,B:1.0e20);")]
  #[case::plain_comments(    NewickDialect::CLASSIC, "([a]A[b]:[c]1[d],[e]:2)[r]Y[f];")]
  #[case::ampersand_classic( NewickDialect::CLASSIC, "(A[&a=1],B);")]
  #[case::beast(             NewickDialect::BEAST,   "[&R]((A[&rate=1.5,s=\"x,y\"]:[&r=1]1[c],B)[&posterior=0.9]:1,C);")]
  #[case::beast_arrays(      NewickDialect::BEAST,   "(A[&a={1,{2,\"x\"},{}},c=#ff0080,b=TRUE],B);")]
  #[case::mrbayes(           NewickDialect::MRBAYES, "[&U](A:1[&B TK02Brlens 0.1],B[&E ibr 2: 0.1]);")]
  #[case::nhx(               NewickDialect::NHX,     "(A[&&NHX:S=human:T=9606:C=1.2.3:D:Ev=1>2]:1,B);")]
  #[case::enewick(           NewickDialect::ENEWICK, "((A)x##LGT1,(x#LGT1,B));")]
  #[case::enewick_unnamed(   NewickDialect::ENEWICK, "((A)#H1:1,(#H1:2,B));")]
  #[case::rich(              NewickDialect::RICH,    "[&U][&W 0.5]((A:1:90,(B)#H1:::0.3),(#H1:::0.7,C)):0.1;")]
  #[case::rich_field_comment(NewickDialect::RICH,    "(A:1[c]:[d]80,B);")]
  #[trace]
  fn test_write_reproduces_canonical_text(#[case] dialect: NewickDialect, #[case] text: &str) {
    assert_eq!(text, rewrite(text, dialect, &NewickWriteOptions::new(dialect)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::spaces(              "a b",       false, "'a b'")]
  #[case::parentheses(         "node (1)",  false, "'node (1)'")]
  #[case::apostrophe(          "it's",      false, "'it''s'")]
  #[case::simple(              "simple_n",  false, "simple_n")]
  #[case::hybrid_shaped(       "A#1",       false, "'A#1'")]
  #[case::hash_without_digits( "A#x",       false, "A#x")]
  #[case::empty(               "",          false, "''")]
  #[case::numeric_leaf(        "123",       false, "123")]
  #[case::numeric_internal(    "123",       true,  "'123'")]
  #[case::support_shaped(      "80/95",     true,  "'80/95'")]
  #[case::infinity_internal(   "inf",       true,  "inf")]
  #[trace]
  fn test_write_quotes_names_by_grammar(#[case] name: &str, #[case] internal: bool, #[case] expected: &str) {
    let graph = if internal {
      star(NewickNodeData::new().with_name(name), vec![NewickNodeData::new()])
    } else {
      NewickGraph::new(NewickNodeData::new().with_name(name))
    };

    let written = newick_to_string(&graph, &NewickWriteOptions::default()).unwrap();

    let expected_text = if internal { format!("(){expected};") } else { format!("{expected};") };
    assert_eq!((expected_text, true), (written.clone(), graph.eq_ordered(&reread(&written, NewickDialect::CLASSIC))));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::always(       NewickWriteOptions { quoting: Quoting::Always, ..NewickWriteOptions::default() },   "('a b','c');")]
  #[case::underscores(  NewickWriteOptions { spaces: Spaces::Underscore, ..NewickWriteOptions::default() }, "(a_b,c);")]
  #[trace]
  fn test_write_quoting_options(#[case] options: NewickWriteOptions, #[case] expected: &str) {
    let graph = star(NewickNodeData::new(), vec![NewickNodeData::new().with_name("a b"), NewickNodeData::new().with_name("c")]);

    assert_eq!(expected, newick_to_string(&graph, &options).unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::shortest(        NumberFormat::default(),                                                              "(A:0.123456789,B:1);")]
  #[case::significant(     NumberFormat { significant_digits: Some(3), ..NumberFormat::default() },              "(A:0.123,B:1);")]
  #[case::decimal(         NumberFormat { decimal_digits: Some(2), ..NumberFormat::default() },                  "(A:0.12,B:1);")]
  #[case::point_zero(      NumberFormat { point_zero: true, ..NumberFormat::default() },                         "(A:0.123456789,B:1.0);")]
  #[case::rounded_zeros(   NumberFormat { significant_digits: Some(3), ..NumberFormat::default() },              "(A:0.123,B:1);")]
  #[trace]
  fn test_write_number_format(#[case] numbers: NumberFormat, #[case] expected: &str) {
    let options = NewickWriteOptions {
      numbers,
      ..NewickWriteOptions::default()
    };

    assert_eq!(expected, rewrite("(A:0.123456789,B:1);", NewickDialect::CLASSIC, &options));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::positional(   0.10004,  Some(3), false, "0.1")]
  #[case::positional_pz(0.10004,  Some(3), true,  "0.1")]
  #[case::integer_pz(   2.0004,   Some(3), true,  "2.0")]
  #[case::exponent(     6.0004e-5, Some(3), false, "6.00e-5")]
  #[case::kept_digits(  0.1234,   Some(3), false, "0.123")]
  #[trace]
  fn test_write_number_format_trims_rounded_zeros(#[case] value: f64, #[case] digits: Option<u8>, #[case] point_zero: bool, #[case] expected: &str) {
    let format = NumberFormat {
      significant_digits: digits,
      decimal_digits: None,
      point_zero,
    };

    assert_eq!(expected, format.format(value).unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::small_negative(  -0.00001,  "(:-0);")]
  #[case::small_positive(  0.00049,   "(:0);")]
  #[case::rounds_up(       0.0006,    "(:0.001);")]
  #[trace]
  fn test_write_decimal_digits_round_small_numbers(#[case] length: f64, #[case] expected: &str) {
    let graph = star_with_edge(NewickEdgeData::new().with_length(length));
    let options = NewickWriteOptions {
      numbers: NumberFormat { decimal_digits: Some(3), ..NumberFormat::default() },
      ..NewickWriteOptions::default()
    };

    assert_eq!(expected, newick_to_string(&graph, &options).unwrap());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::label_to_field(    "((A,B)95,C);",          NewickDialect::RICH,    SupportPlacement::Field,                           "((A,B)::95,C);")]
  #[case::field_to_label(    "((A,B):1:95,C);",       NewickDialect::RICH,    SupportPlacement::Label,                           "((A,B)95:1,C);")]
  #[case::field_in_classic(  "((A,B):1:95,C);",       NewickDialect::CLASSIC, SupportPlacement::Source,                          "((A,B):1,C);")]
  #[case::annotation_beast(  "((A,B)0.9,C);",         NewickDialect::BEAST,   SupportPlacement::Annotation("posterior".to_owned()), "((A,B)[&posterior=0.9],C);")]
  #[case::annotation_nhx(    "((A,B)90/80,C);",       NewickDialect::NHX,     SupportPlacement::Annotation("s".to_owned()),      "((A,B)[&&NHX:s=90>80],C);")]
  #[trace]
  fn test_write_support_placement(#[case] input: &str, #[case] dialect: NewickDialect, #[case] support: SupportPlacement, #[case] expected: &str) {
    let options = NewickWriteOptions {
      support,
      ..NewickWriteOptions::new(dialect)
    };

    assert_eq!(expected, rewrite(input, NewickDialect::RICH, &options));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::leaf_label(       star_with_edge(NewickEdgeData::new().with_support(vec![0.9], SupportSource::Label)),  NewickWriteOptions::default(),  "When writing Newick: When writing the branch above node 1: A leaf cannot carry support in its label, because a leaf label is read as a name")]
  #[case::named_label(      named_internal_with_support(),                                                       NewickWriteOptions::default(),  "When writing Newick: When writing the branch above node 1 ('X'): The label of the node holds its name, so it cannot also hold the support of the branch above")]
  #[case::field_classic(    star_with_edge(NewickEdgeData::new()),                                               NewickWriteOptions { support: SupportPlacement::Field, ..NewickWriteOptions::default() }, "When writing Newick: Support in a colon field needs the rich structure, not classic,plain")]
  #[case::annotation_rich(  star_with_edge(NewickEdgeData::new()),                                               NewickWriteOptions { support: SupportPlacement::Annotation("p".to_owned()), ..NewickWriteOptions::new(NewickDialect::RICH) }, "When writing Newick: Support in an annotation needs beast or nhx annotations, not rich,plain")]
  #[case::two_in_field(     star_with_edge(NewickEdgeData::new().with_support(vec![1.0, 2.0], SupportSource::Field)), NewickWriteOptions::new(NewickDialect::RICH), "When writing Newick: When writing the branch above node 1: A colon field holds one support value, but the branch has 2")]
  #[case::infinite_length(  star_with_edge(NewickEdgeData::new().with_length(f64::INFINITY)),                    NewickWriteOptions::default(),  "When writing Newick: When writing the branch above node 1: Newick cannot represent the number inf")]
  #[case::nan_support(      named_root_support(f64::NAN),                                                        NewickWriteOptions::default(),  "When writing Newick: When writing node 0: Newick cannot represent the number NaN")]
  #[case::zero_digits(      star_with_edge(NewickEdgeData::new().with_length(1.0)),                              NewickWriteOptions { numbers: NumberFormat { significant_digits: Some(0), ..NumberFormat::default() }, ..NewickWriteOptions::default() }, "When writing Newick: When writing the branch above node 1: The number of significant digits must be at least 1")]
  #[trace]
  fn test_write_rejects_unrepresentable(#[case] graph: NewickGraph, #[case] options: NewickWriteOptions, #[case] expected: &str) {
    assert_eq!(expected, write_error(&graph, &options));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::before(   BranchAnnotations::BeforeLength, "(A:[c][&p=1]1,B:[&q=2]);")]
  #[case::after(    BranchAnnotations::AfterLength,  "(A:[c]1[&p=1],B:[&q=2]);")]
  #[case::recorded( BranchAnnotations::Recorded,     "(A:[c]1[&p=1],B:[&q=2]);")]
  #[trace]
  fn test_write_branch_annotation_position(#[case] branch_annotations: BranchAnnotations, #[case] expected: &str) {
    let options = NewickWriteOptions {
      branch_annotations,
      ..NewickWriteOptions::new(NewickDialect::BEAST)
    };

    assert_eq!(expected, rewrite("(A:[c]1[&p=1],B:[&q=2]);", NewickDialect::BEAST, &options));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::string(       vec![("s", NewickValue::String("usa".to_owned()))],              "[&s=\"usa\"]")]
  #[case::quote(        vec![("s", NewickValue::String("a\"b".to_owned()))],             "[&s=\"a\"\"b\"]")]
  #[case::number_text(  vec![("d", NewickValue::NumberText("2020.50".to_owned()))],      "[&d=2020.50]")]
  #[case::huge_number(  vec![("k", NewickValue::Number(1e300))],                         "[&k=1.0e300]")]
  #[case::quoted_key(   vec![("posterior prob", NewickValue::Number(0.5))],              "[&\"posterior prob\"=0.5]")]
  #[case::true_key(     vec![("TRUE", NewickValue::Boolean(false))],                     "[&TRUE=FALSE]")]
  #[case::empty_key(    vec![("", NewickValue::Boolean(true))],                          "[&\"\"=TRUE]")]
  #[case::caller_order( vec![("m", NewickValue::String("A55G".to_owned())), ("d", NewickValue::NumberText("2003.80".to_owned()))], "[&m=\"A55G\",d=2003.80]")]
  #[trace]
  fn test_write_beast_values(#[case] pairs: Vec<(&str, NewickValue)>, #[case] expected: &str) {
    let graph = commented_leaf(NewickComment::Beast(helpers::pairs(pairs)));

    let written = newick_to_string(&graph, &NewickWriteOptions::new(NewickDialect::BEAST)).unwrap();

    assert_eq!((format!("A{expected};"), true), (written.clone(), graph.eq_ordered(&reread(&written, NewickDialect::BEAST))));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::reserved_colon(  vec![("k", NewickValue::String("a:b".to_owned()))],   "When writing Newick: When writing node 0 ('A'): NHX cannot hold the value \"a:b\" of the tag k: a value contains none of ':', '=', '[', ']', '>'")]
  #[case::reserved_key(    vec![("a=b", NewickValue::String("x".to_owned()))],   "When writing Newick: When writing node 0 ('A'): NHX cannot hold the tag \"a=b\": a tag is not empty and contains none of ':', '=', '[', ']', '>'")]
  #[case::false_value(     vec![("k", NewickValue::Boolean(false))],             "When writing Newick: When writing node 0 ('A'): The NHX tag k needs text, and NHX cannot hold the value Boolean(false) there")]
  #[case::typed_mismatch(  vec![("B", NewickValue::String("high".to_owned()))],  "When writing Newick: When writing node 0 ('A'): The NHX tag B needs a decimal number, and NHX cannot hold the value String(\"high\") there")]
  #[case::integer_tag(     vec![("T", NewickValue::Number(1.5))],                "When writing Newick: When writing node 0 ('A'): The NHX tag T needs an integer, but the value is 1.5")]
  #[case::single_array(    vec![("k", NewickValue::Array(vec![NewickValue::String("x".to_owned())].into()))], "When writing Newick: When writing node 0 ('A'): The NHX tag k needs text, and NHX cannot hold the value Array(NewickArray([String(\"x\")])) there")]
  #[trace]
  fn test_write_nhx_rejects(#[case] pairs: Vec<(&str, NewickValue)>, #[case] expected: &str) {
    let graph = commented_leaf(NewickComment::Nhx(helpers::pairs(pairs)));

    assert_eq!(expected, write_error(&graph, &NewickWriteOptions::new(NewickDialect::NHX)));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unbalanced(       NewickDialect::CLASSIC, "a]b",  "When writing Newick: When writing node 0 ('A'): The comment text \"a]b\" cannot be written as a comment, because its brackets do not balance")]
  #[case::ampersand_beast(  NewickDialect::BEAST,   "&a=1", "When writing Newick: When writing node 0 ('A'): The comment Plain(\"&a=1\") cannot be written in the classic,beast dialect, which would read it as an annotation")]
  #[case::ampersand_rich(   NewickDialect::new(NewickStructure::Rich, NewickAnnotations::Beast), "&R", "When writing Newick: When writing node 0 ('A'): The comment Plain(\"&R\") cannot be written in the rich,beast dialect, which would read it as an annotation")]
  #[trace]
  fn test_write_rejects_plain_comment(#[case] dialect: NewickDialect, #[case] text: &str, #[case] expected: &str) {
    let graph = commented_leaf(NewickComment::Plain(text.to_owned()));

    assert_eq!(expected, write_error(&graph, &NewickWriteOptions::new(dialect)));
  }

  #[test]
  fn test_write_hybrid_name_ending_in_hash_is_quoted() {
    let text = "(A,(B)'x#'#H1,('x#'#H1,C));";

    assert_eq!(
      text,
      rewrite(
        text,
        NewickDialect::ENEWICK,
        &NewickWriteOptions::new(NewickDialect::ENEWICK)
      )
    );
  }

  #[test]
  fn test_write_indent() {
    let options = NewickWriteOptions {
      indent: Some(2),
      ..NewickWriteOptions::default()
    };

    let written = rewrite("((A:1,B:2)C:3,D)E;", NewickDialect::CLASSIC, &options);

    assert_eq!(
      indoc::indoc! {"
        (
          (
            A:1,
            B:2
          )C:3,
          D
        )E;"},
      written
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::kept(     true,  "(A,B)95:0.5;")]
  #[case::dropped(  false, "(A,B);")]
  #[trace]
  fn test_write_root_edge_option(#[case] root_edge: bool, #[case] expected: &str) {
    let options = NewickWriteOptions {
      root_edge,
      ..NewickWriteOptions::default()
    };

    assert_eq!(expected, rewrite("(A,B)95:0.5;", NewickDialect::CLASSIC, &options));
  }

  #[test]
  fn test_write_several_trees_one_per_line() {
    let first = reread("(A,B);", NewickDialect::CLASSIC);
    let second = reread("(C);", NewickDialect::CLASSIC);
    let mut buffer = Vec::new();

    write_newick_trees(&mut buffer, [&first, &second], &NewickWriteOptions::default()).unwrap();

    assert_eq!("(A,B);\n(C);\n", String::from_utf8(buffer).unwrap());
  }

  #[test]
  fn test_write_deep_tree_on_small_stack() {
    let input = caterpillar(100_000);
    let expected = input.clone();

    let written = on_small_stack(move || {
      let graph = read_with(&input, &NewickReadOptions::default()).graph;
      newick_to_string(&graph, &NewickWriteOptions::default()).unwrap()
    });

    assert!(written == expected, "the written deep tree differs from the input");
  }

  #[test]
  fn test_write_deep_array_on_small_stack() {
    let input = format!("A[&a={}1{}];", "{".repeat(100_000), "}".repeat(100_000));
    let expected = input.clone();

    let written = on_small_stack(move || {
      rewrite(
        &input,
        NewickDialect::BEAST,
        &NewickWriteOptions::new(NewickDialect::BEAST),
      )
    });

    assert!(written == expected, "the written deep array differs from the input");
  }

  mod helpers {
    use crate::__tests__::test_read_basic::tests::helpers::read_with;
    use crate::dialect::NewickDialect;
    use crate::model::comment::{LabelSide, NewickComment, NodeComment};
    use crate::model::data::{NewickEdgeData, NewickNodeData, SupportSource};
    use crate::model::graph::NewickGraph;
    use crate::model::value::NewickValue;
    use crate::read::options::NewickReadOptions;
    use crate::write::newick::newick_to_string;
    use crate::write::options::NewickWriteOptions;

    pub(super) fn reread(text: &str, dialect: NewickDialect) -> NewickGraph {
      let options = NewickReadOptions {
        dialect,
        ..NewickReadOptions::default()
      };
      read_with(text, &options).graph
    }

    pub(super) fn rewrite(text: &str, dialect: NewickDialect, options: &NewickWriteOptions) -> String {
      newick_to_string(&reread(text, dialect), options).unwrap()
    }

    pub(super) fn write_error(graph: &NewickGraph, options: &NewickWriteOptions) -> String {
      format!("{:#}", newick_to_string(graph, options).unwrap_err())
    }

    pub(super) fn star(root: NewickNodeData, leaves: Vec<NewickNodeData>) -> NewickGraph {
      let mut graph = NewickGraph::new(root);
      for leaf in leaves {
        graph.add_child(0, NewickEdgeData::new(), leaf).unwrap();
      }
      graph
    }

    pub(super) fn pairs(pairs: Vec<(&str, NewickValue)>) -> Vec<(String, NewickValue)> {
      pairs.into_iter().map(|(key, value)| (key.to_owned(), value)).collect()
    }

    pub(super) fn star_with_edge(edge: NewickEdgeData) -> NewickGraph {
      let mut graph = NewickGraph::new(NewickNodeData::new());
      graph.add_child(0, edge, NewickNodeData::new()).unwrap();
      graph
    }

    pub(super) fn named_internal_with_support() -> NewickGraph {
      let mut graph = NewickGraph::new(NewickNodeData::new());
      let inner = graph
        .add_child(
          0,
          NewickEdgeData::new().with_support(vec![0.9], SupportSource::Label),
          NewickNodeData::new().with_name("X"),
        )
        .unwrap();
      graph
        .add_child(inner, NewickEdgeData::new(), NewickNodeData::new())
        .unwrap();
      graph
    }

    pub(super) fn named_root_support(value: f64) -> NewickGraph {
      let mut graph = star_with_edge(NewickEdgeData::new());
      graph.root_edge_mut().set_support(vec![value], SupportSource::Label);
      graph
    }

    pub(super) fn commented_leaf(comment: NewickComment) -> NewickGraph {
      NewickGraph::new(
        NewickNodeData::new()
          .with_name("A")
          .with_comment(NodeComment::new(LabelSide::AfterLabel, comment)),
      )
    }
  }
}
