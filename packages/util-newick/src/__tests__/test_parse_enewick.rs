#[cfg(test)]
mod tests {
  use crate::types::{NewickHybrid, NewickLabel};
  use helpers::{hybrid_nodes, parse_enewick, parse_enewick_error};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[test]
  fn test_parse_enewick_quoted_label_is_never_hybrid() {
    let g = parse_enewick("('A#1',B);").unwrap();

    assert_eq!((Some("A#1"), 0), (g.nodes[0].name(), hybrid_nodes(&g).len()));
  }

  #[test]
  fn test_parse_enewick_quoted_name_with_hybrid_tag() {
    let g = parse_enewick("((C)'x y'#H1,('x y'#H1,D));").unwrap();

    let expected = vec![(
      Some(NewickLabel::Name("x y".to_owned())),
      NewickHybrid {
        kind: Some("H".to_owned()),
        index: 1,
      },
    )];
    assert_eq!(expected, hybrid_nodes(&g));
  }

  #[test]
  fn test_parse_enewick_non_ascii_digits_are_not_an_index() {
    let g = parse_enewick("(A#\u{0663},B);").unwrap();

    assert_eq!((Some("A#\u{0663}"), 0), (g.nodes[0].name(), hybrid_nodes(&g).len()));
  }

  #[test]
  fn test_parse_enewick_label_of_later_occurrence_is_kept() {
    let g = parse_enewick("((C)#H1,(X#H1,D));").unwrap();

    let expected = vec![(
      Some(NewickLabel::Name("X".to_owned())),
      NewickHybrid {
        kind: Some("H".to_owned()),
        index: 1,
      },
    )];
    assert_eq!(expected, hybrid_nodes(&g));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::different_labels(  "((C)x#H1,(y#H1,D));",   "The occurrences of the hybrid node #H1 have different labels, at node 1 ('x')")]
  #[case::children_twice(    "((A)#H1,(B)#H1);",      "The hybrid node #H1 has children in more than one of its occurrences")]
  #[case::parallel_edges(    "((C)#H1,#H1);",         "node 2 has more than one edge to node 1")]
  #[case::self_loop(         "(A,#H1)#H1;",           "The root, node 1, has a parent")]
  #[case::cycle(             "((#H2)#H1,(#H1)#H2);",  "The graph contains a cycle")]
  #[case::index_too_large(   "(A#H4294967296,B);",    "When parsing the hybrid node index in 'A#H4294967296': number too large to fit in target type")]
  #[trace]
  fn test_parse_enewick_invalid_network_is_error(#[case] input: &str, #[case] expected: &str) {
    assert_eq!(format!("Failed to parse Newick string: {expected}"), parse_enewick_error(input));
  }

  mod helpers {
    use crate::parse::newick_from_string;
    use crate::types::{NewickGraph, NewickHybrid, NewickLabel, NewickReadOptions};

    pub(super) fn parse_enewick(input: &str) -> eyre::Result<NewickGraph> {
      newick_from_string(input, &NewickReadOptions { enewick: true })
    }

    pub(super) fn parse_enewick_error(input: &str) -> String {
      format!("{:#}", parse_enewick(input).unwrap_err())
    }

    pub(super) fn hybrid_nodes(graph: &NewickGraph) -> Vec<(Option<NewickLabel>, NewickHybrid)> {
      graph
        .nodes
        .iter()
        .filter_map(|node| Some((node.label.clone(), node.hybrid.clone()?)))
        .collect()
    }
  }
}
