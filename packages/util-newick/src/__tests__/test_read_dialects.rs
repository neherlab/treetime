#[cfg(test)]
mod tests {
  use crate::dialect::{NewickAnnotations, NewickDialect, NewickStructure};
  use crate::read::error::NewickErrorKind;
  use crate::read::options::{NewickReadOptions, ReadMode};
  use crate::read::stream::newick_from_str;
  use helpers::{Shape, shape};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  const COALRE: &str = "((A[&segments={0,1}]:1,(B[&segments={0,1}]:1)#H0[&segments={0}]:0.5)[&segments={0,1}]:1,(#H0[&segments={1}]:0.7,C[&segments={0,1}]:1)[&segments={0,1}]:1);";
  const RICH_BEAST: &str = "[&U]((A[&x=1]:1:90,(B)#H1[&gamma=0.3]:::0.3),(#H1[&gamma=0.7]:::0.7,C));";
  const GISAID: &str = "(EPI_ISL#402124:1,B:1);";
  const ROOTING: &str = "[&U]((A,B),C);";
  const NHX: &str = "((A[&&NHX:S=human]:1,B:1):1,C:2);";

  const RICH_BEAST_DIALECT: NewickDialect = NewickDialect::new(NewickStructure::Rich, NewickAnnotations::Beast);
  const RICH_PLAIN_DIALECT: NewickDialect = NewickDialect::new(NewickStructure::Rich, NewickAnnotations::Plain);

  #[rustfmt::skip]
  #[rstest]
  #[case::coalre_enewick_beast(  COALRE,     NewickDialect::ENEWICK_BEAST, Shape {
    rooted: None,
    hybrid_parents: vec![2],
    labels: vec!["-", "-", "-", "-[#H0]", "A", "B", "C"],
    node_comments: vec!["beast:segments={0,1}"; 5],
    occurrence_comments: vec![":0.5 beast:segments={0}", ":0.7 beast:segments={1}"],
  })]
  #[case::coalre_beast(          COALRE,     NewickDialect::BEAST,         Shape {
    rooted: None,
    hybrid_parents: vec![],
    labels: vec!["#H0", "#H0", "-", "-", "-", "A", "B", "C"],
    node_comments: vec!["beast:segments={0,1}", "beast:segments={0,1}", "beast:segments={0,1}", "beast:segments={0,1}", "beast:segments={0,1}", "beast:segments={0}", "beast:segments={1}"],
    occurrence_comments: vec![],
  })]
  #[case::coalre_enewick(        COALRE,     NewickDialect::ENEWICK,       Shape {
    rooted: None,
    hybrid_parents: vec![2],
    labels: vec!["-", "-", "-", "-[#H0]", "A", "B", "C"],
    node_comments: vec!["plain:&segments={0,1}"; 5],
    occurrence_comments: vec![":0.5 plain:&segments={0}", ":0.7 plain:&segments={1}"],
  })]
  #[case::rich_beast(            RICH_BEAST, RICH_BEAST_DIALECT,           Shape {
    rooted: Some(false),
    hybrid_parents: vec![2],
    labels: vec!["-", "-", "-", "-[#H1]", "A", "B", "C"],
    node_comments: vec!["beast:x=1"],
    occurrence_comments: vec!["p=0.3 beast:gamma=0.3", "p=0.7 beast:gamma=0.7"],
  })]
  #[case::gisaid_beast(          GISAID,     NewickDialect::BEAST,         Shape {
    rooted: None,
    hybrid_parents: vec![],
    labels: vec!["-", "B", "EPI_ISL#402124"],
    node_comments: vec![],
    occurrence_comments: vec![],
  })]
  #[case::gisaid_enewick(        GISAID,     NewickDialect::ENEWICK,       Shape {
    rooted: None,
    hybrid_parents: vec![1],
    labels: vec!["-", "B", "EPI_ISL[#402124]"],
    node_comments: vec![],
    occurrence_comments: vec![],
  })]
  #[case::rooting_rich(          ROOTING,    NewickDialect::RICH,          Shape {
    rooted: Some(false),
    hybrid_parents: vec![],
    labels: vec!["-", "-", "A", "B", "C"],
    node_comments: vec![],
    occurrence_comments: vec![],
  })]
  #[case::rooting_classic(       ROOTING,    NewickDialect::CLASSIC,       Shape {
    rooted: None,
    hybrid_parents: vec![],
    labels: vec!["-", "-", "A", "B", "C"],
    node_comments: vec!["plain:&U"],
    occurrence_comments: vec![],
  })]
  #[trace]
  fn test_read_dialects_pair(#[case] input: &str, #[case] dialect: NewickDialect, #[case] expected: Shape) {
    let options = NewickReadOptions {
      dialect,
      ..NewickReadOptions::default()
    };

    let tree = newick_from_str(input, &options).unwrap();

    assert_eq!(expected, shape(&tree.graph));
  }

  #[test]
  fn test_read_dialects_rich_beast_support() {
    let options = NewickReadOptions {
      dialect: RICH_BEAST_DIALECT,
      ..NewickReadOptions::default()
    };
    let tree = newick_from_str(RICH_BEAST, &options).unwrap();

    let support: Vec<Vec<f64>> = tree
      .graph
      .edges()
      .map(|(_, edge)| edge.data().support().to_vec())
      .filter(|support| !support.is_empty())
      .collect();

    assert_eq!(vec![vec![90.0]], support);
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rich_plain(  RICH_PLAIN_DIALECT, ROOTING)]
  #[case::rich_beast(  RICH_BEAST_DIALECT, ROOTING)]
  #[trace]
  fn test_read_dialects_rich_reads_rooting(#[case] dialect: NewickDialect, #[case] input: &str) {
    let options = NewickReadOptions {
      dialect,
      ..NewickReadOptions::default()
    };

    let tree = newick_from_str(input, &options).unwrap();

    assert_eq!(Some(false), tree.graph.rooted());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::strict(    ReadMode::Strict,   Err(NewickErrorKind::Annotation))]
  #[case::tolerant(  ReadMode::Tolerant, Ok(1))]
  #[trace]
  fn test_read_dialects_nhx_as_beast(#[case] mode: ReadMode, #[case] expected: Result<usize, NewickErrorKind>) {
    let options = NewickReadOptions {
      dialect: NewickDialect::BEAST,
      mode,
      ..NewickReadOptions::default()
    };

    let actual = newick_from_str(NHX, &options).map(|tree| tree.warnings.len()).map_err(|error| error.kind);

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_read_dialects_every_pair_reads_a_plain_tree() {
    let read: Vec<String> = NewickDialect::pairs()
      .filter(|&dialect| {
        let options = NewickReadOptions {
          dialect,
          ..NewickReadOptions::default()
        };
        newick_from_str("((A:1,B:2)C:3,D:4);", &options).is_ok()
      })
      .map(|dialect| dialect.to_string())
      .collect();

    let expected: Vec<String> = NewickDialect::pairs().map(|dialect| dialect.to_string()).collect();
    assert_eq!(expected, read);
  }

  mod helpers {
    use crate::__tests__::test_read_basic::tests::helpers::edge_summary;
    use crate::model::comment::NewickComment;
    use crate::model::graph::NewickGraph;
    use crate::model::value::NewickValue;

    #[derive(Debug, PartialEq)]
    pub(super) struct Shape {
      pub(super) rooted: Option<bool>,
      pub(super) hybrid_parents: Vec<usize>,
      pub(super) labels: Vec<&'static str>,
      pub(super) node_comments: Vec<&'static str>,
      pub(super) occurrence_comments: Vec<&'static str>,
    }

    pub(super) fn shape(graph: &NewickGraph) -> ActualShape {
      let mut labels = Vec::new();
      let mut node_comments = Vec::new();
      let mut hybrid_parents = Vec::new();
      for (node, data) in graph.nodes() {
        let name = data.name().unwrap_or("-");
        labels.push(match data.hybrid() {
          Some(hybrid) => format!("{name}[{}]", hybrid.tag(false)),
          None => name.to_owned(),
        });
        if data.hybrid().is_some() {
          hybrid_parents.push(graph.parent_edges(node).len());
        }
        node_comments.extend(data.comments().iter().map(|comment| render(&comment.comment)));
      }
      let occurrence_comments = graph
        .edges()
        .map(|(_, edge)| edge.data())
        .chain([graph.root_edge()])
        .flat_map(|edge| {
          edge
            .occurrence_comments()
            .iter()
            .map(move |comment| format!("{} {}", edge_summary(edge), render(&comment.comment)))
        })
        .collect();
      ActualShape {
        rooted: graph.rooted(),
        hybrid_parents,
        labels: sorted(labels),
        node_comments: sorted(node_comments),
        occurrence_comments: sorted(occurrence_comments),
      }
    }

    #[derive(Debug)]
    pub(super) struct ActualShape {
      rooted: Option<bool>,
      hybrid_parents: Vec<usize>,
      labels: Vec<String>,
      node_comments: Vec<String>,
      occurrence_comments: Vec<String>,
    }

    impl PartialEq<ActualShape> for Shape {
      fn eq(&self, other: &ActualShape) -> bool {
        self.rooted == other.rooted
          && self.hybrid_parents == other.hybrid_parents
          && self.labels == other.labels
          && self.node_comments == other.node_comments
          && self.occurrence_comments == other.occurrence_comments
      }
    }

    fn sorted(mut items: Vec<String>) -> Vec<String> {
      items.sort();
      items
    }

    fn render(comment: &NewickComment) -> String {
      match comment {
        NewickComment::Plain(text) => format!("plain:{text}"),
        NewickComment::Beast(pairs) => {
          let pairs: Vec<String> = pairs
            .iter()
            .map(|(key, value)| format!("{key}={}", value_text(value)))
            .collect();
          format!("beast:{}", pairs.join(","))
        },
        NewickComment::Nhx(_) | NewickComment::MrBayesMcmc(_) => format!("{comment:?}"),
      }
    }

    fn value_text(value: &NewickValue) -> String {
      match value {
        NewickValue::Number(number) => number.to_string(),
        NewickValue::Array(values) => {
          let items: Vec<String> = values.as_slice().iter().map(value_text).collect();
          format!("{{{}}}", items.join(","))
        },
        NewickValue::NumberText(_) | NewickValue::String(_) | NewickValue::Boolean(_) | NewickValue::Color(_) => {
          format!("{value:?}")
        },
      }
    }
  }
}
