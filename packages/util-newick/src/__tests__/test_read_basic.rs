#[cfg(test)]
pub(crate) mod tests {
  use crate::dialect::NewickDialect;
  use crate::read::options::{InternalLabel, NewickReadOptions};
  use helpers::{read, read_with, summary};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::no_names(              "(,,(,));",                            vec!["- ", "- ", "- ", "- ", "- ", "- "])]
  #[case::leaf_names(            "(A,B,(C,D));",                        vec!["- ", "A ", "B ", "- ", "C ", "D "])]
  #[case::all_names(             "(A,B,(C,D)E)F;",                      vec!["F ", "A ", "B ", "E ", "C ", "D "])]
  #[case::distances_leaf_names(  "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);",    vec!["- ", "A :0.1", "B :0.2", "- :0.5", "C :0.3", "D :0.4"])]
  #[case::distances_all_names(   "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;",  vec!["F ", "A :0.1", "B :0.2", "E :0.5", "C :0.3", "D :0.4"])]
  #[case::single_leaf(           "A;",                                  vec!["A "])]
  #[case::lone_semicolon(        ";",                                   vec!["- "])]
  #[case::single_child(          "((A));",                              vec!["- ", "- ", "A "])]
  #[trace]
  fn test_read_basic_felsenstein_examples(#[case] input: &str, #[case] expected: Vec<&str>) {
    assert_eq!(expected, summary(&read(input).graph));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::scientific(       "(A:1.5e-3,B:2E4);",         vec!["- ", "A :0.0015", "B :20000"])]
  #[case::negative(         "(A:-0.01,B:0.2);",          vec!["- ", "A :-0.01", "B :0.2"])]
  #[case::leading_dot(      "(A:.5,B:2);",               vec!["- ", "A :0.5", "B :2"])]
  #[case::plus_sign(        "(A:+1,B:2);",               vec!["- ", "A :1", "B :2"])]
  #[case::overflow(         "(A:1e400,B:1);",            vec!["- ", "A :inf", "B :1"])]
  #[case::root_length(      "(A:0.1,B:0.2):0.5;",        vec!["- :0.5", "A :0.1", "B :0.2"])]
  #[case::empty_field(      "(A:,B);",                   vec!["- ", "A ", "B "])]
  #[case::whitespace(       "  ( A : 0.1 , B : 0.2 ) ; ", vec!["- ", "A :0.1", "B :0.2"])]
  #[case::newlines(         "(\n  A:0.1,\n  B:0.2\n);\n", vec!["- ", "A :0.1", "B :0.2"])]
  #[trace]
  fn test_read_basic_branch_lengths(#[case] input: &str, #[case] expected: Vec<&str>) {
    assert_eq!(expected, summary(&read(input).graph));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::quoted(              "'node with spaces';",   "node with spaces")]
  #[case::escaped_quote(       "'it''s a name';",       "it's a name")]
  #[case::empty_quoted(        "'';",                   "")]
  #[case::apostrophe_inside(   "it's;",                 "it's")]
  #[case::hash(                "A#1;",                  "A#1")]
  #[case::underscores_kept(    "Homo_sapiens;",         "Homo_sapiens")]
  #[case::number_leaf(         "0.999;",                "0.999")]
  #[case::quoted_brackets(     "'a[b]:c';",             "a[b]:c")]
  #[trace]
  fn test_read_basic_names(#[case] input: &str, #[case] expected: &str) {
    let tree = read(input);

    assert_eq!(Some(expected), tree.graph.node(tree.graph.root()).name());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::unquoted(  "(Homo_sapiens,'Pan_troglodytes');", vec!["- ", "Homo sapiens ", "Pan_troglodytes "])]
  #[trace]
  fn test_read_basic_underscores_as_spaces(#[case] input: &str, #[case] expected: Vec<&str>) {
    let options = NewickReadOptions {
      underscores_as_spaces: true,
      ..NewickReadOptions::default()
    };

    assert_eq!(expected, summary(&read_with(input, &options).graph));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::auto_float(        InternalLabel::Auto,    "((A,B)0.999:0.3,C);",  vec!["- ", "- :0.3 support=0.999(label)", "A ", "B ", "C "])]
  #[case::auto_integer(      InternalLabel::Auto,    "((A,B)100,C);",        vec!["- ", "- support=100(label)", "A ", "B ", "C "])]
  #[case::auto_multi(        InternalLabel::Auto,    "((A,B)80.5/95,C);",    vec!["- ", "- support=80.5/95(label)", "A ", "B ", "C "])]
  #[case::auto_root(         InternalLabel::Auto,    "(A,B)0.5;",            vec!["- support=0.5(label)", "A ", "B "])]
  #[case::auto_name(         InternalLabel::Auto,    "((A,B)clade,C);",      vec!["- ", "clade ", "A ", "B ", "C "])]
  #[case::auto_quoted(       InternalLabel::Auto,    "((A,B)'95',C);",       vec!["- ", "95 ", "A ", "B ", "C "])]
  #[case::auto_inf_is_name(  InternalLabel::Auto,    "((A,B)inf,C);",        vec!["- ", "inf ", "A ", "B ", "C "])]
  #[case::auto_leaf_number(  InternalLabel::Auto,    "(0.999:0.1,B);",       vec!["- ", "0.999 :0.1", "B "])]
  #[case::name_mode(         InternalLabel::Name,    "((A,B)95,C);",         vec!["- ", "95 ", "A ", "B ", "C "])]
  #[case::support_mode(      InternalLabel::Support, "((A,B)95,C);",         vec!["- ", "- support=95(label)", "A ", "B ", "C "])]
  #[trace]
  fn test_read_basic_internal_labels(#[case] mode: InternalLabel, #[case] input: &str, #[case] expected: Vec<&str>) {
    let options = NewickReadOptions {
      internal_label: mode,
      ..NewickReadOptions::default()
    };

    assert_eq!(expected, summary(&read_with(input, &options).graph));
  }

  #[test]
  fn test_read_basic_support_mode_rejects_name() {
    let options = NewickReadOptions {
      internal_label: InternalLabel::Support,
      ..NewickReadOptions::default()
    };

    let actual = crate::read::stream::newick_from_str("((A,B)clade,C);", &options).unwrap_err();

    assert_eq!(
      r#"line 1, column 7: The internal label "clade" is not a support value"#,
      actual.to_string()
    );
  }

  #[test]
  fn test_read_basic_records_dialect() {
    assert_eq!(NewickDialect::Classic, read("(A,B);").dialect);
  }

  pub(crate) mod helpers {
    use crate::model::data::{NewickEdgeData, SupportSource};
    use crate::model::graph::NewickGraph;
    use crate::read::options::{NewickReadOptions, NewickTree};
    use crate::read::stream::newick_from_str;
    use std::fmt::Write;
    use std::thread;

    pub(crate) fn read(input: &str) -> NewickTree {
      read_with(input, &NewickReadOptions::default())
    }

    pub(crate) fn read_with(input: &str, options: &NewickReadOptions) -> NewickTree {
      match newick_from_str(input, options) {
        Ok(tree) => tree,
        Err(error) => panic!("{input:?} does not read: {error}"),
      }
    }

    pub(crate) fn read_error(input: &str, options: &NewickReadOptions) -> String {
      match newick_from_str(input, options) {
        Ok(tree) => panic!("{input:?} reads as {:?}", summary(&tree.graph)),
        Err(error) => error.to_string(),
      }
    }

    pub(crate) fn summary(graph: &NewickGraph) -> Vec<String> {
      graph
        .preorder()
        .map(|node| {
          let data = graph.node(node);
          let tag = data.hybrid().map(|hybrid| hybrid.tag(false)).unwrap_or_default();
          let edges: Vec<String> = if node == graph.root() {
            vec![edge_summary(graph.root_edge())]
          } else {
            graph
              .parent_edges(node)
              .iter()
              .map(|&edge| edge_summary(graph.edge(edge).data()))
              .collect()
          };
          format!("{}{tag} {}", data.name().unwrap_or("-"), edges.join(" | "))
        })
        .collect()
    }

    pub(crate) fn edge_summary(edge: &NewickEdgeData) -> String {
      let mut parts = Vec::new();
      if let Some(length) = edge.branch_length() {
        parts.push(format!(":{length}"));
      }
      if !edge.support().is_empty() {
        let values: Vec<String> = edge.support().iter().map(f64::to_string).collect();
        let source = match edge.support_source() {
          SupportSource::Label => "label",
          SupportSource::Field => "field",
        };
        parts.push(format!("support={}({source})", values.join("/")));
      }
      if let Some(probability) = edge.probability() {
        parts.push(format!("p={probability}"));
      }
      if edge.is_acceptor() {
        parts.push("acceptor".to_owned());
      }
      parts.join(" ")
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
