#[cfg(test)]
mod tests {
  use crate::parse::newick_from_string;
  use crate::types::{NewickReadOptions, NewickWriteOptions, NwkStyle};
  use crate::write::newick_to_string;
  use generators::{AnnotationAlphabet, arb_graph, without_annotations};
  use proptest::prelude::*;

  proptest! {
    #[test]
    fn test_prop_roundtrip_beast(graph in arb_graph(AnnotationAlphabet::Beast)) {
      let written = newick_to_string(&graph, &options(NwkStyle::Beast)).unwrap();
      let parsed = newick_from_string(&written, &NewickReadOptions::default()).unwrap();
      prop_assert!(graph.eq_ordered(&parsed), "Round trip changed the graph.\nWritten: {written}\nBefore: {graph:#?}\nAfter: {parsed:#?}");
    }

    #[test]
    fn test_prop_roundtrip_nhx(graph in arb_graph(AnnotationAlphabet::Nhx)) {
      let written = newick_to_string(&graph, &options(NwkStyle::Nhx)).unwrap();
      let parsed = newick_from_string(&written, &NewickReadOptions::default()).unwrap();
      prop_assert!(graph.eq_ordered(&parsed), "Round trip changed the graph.\nWritten: {written}\nBefore: {graph:#?}\nAfter: {parsed:#?}");
    }

    #[test]
    fn test_prop_roundtrip_plain_drops_only_annotations(graph in arb_graph(AnnotationAlphabet::Beast)) {
      let written = newick_to_string(&graph, &options(NwkStyle::Plain)).unwrap();
      let parsed = newick_from_string(&written, &NewickReadOptions::default()).unwrap();
      let expected = without_annotations(&graph);
      prop_assert!(expected.eq_ordered(&parsed), "Round trip changed the graph.\nWritten: {written}");
    }

    #[test]
    fn test_prop_roundtrip_write_idempotent(graph in arb_graph(AnnotationAlphabet::Beast)) {
      let first = newick_to_string(&graph, &options(NwkStyle::Beast)).unwrap();
      let parsed = newick_from_string(&first, &NewickReadOptions::default()).unwrap();
      let second = newick_to_string(&parsed, &options(NwkStyle::Beast)).unwrap();
      prop_assert_eq!(first, second);
    }
  }

  fn options(style: NwkStyle) -> NewickWriteOptions {
    NewickWriteOptions {
      style,
      significant_digits: None,
      decimal_digits: None,
    }
  }

  mod generators {
    use crate::types::{NewickEdgeData, NewickGraph, NewickLabel, NewickNodeData, NewickValue};
    use proptest::collection::{btree_map, vec};
    use proptest::prelude::*;
    use std::collections::BTreeMap;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum AnnotationAlphabet {
      Beast,
      Nhx,
    }

    #[derive(Clone, Debug)]
    pub(super) struct GenNode {
      node: NewickNodeData,
      children: Vec<(GenNode, NewickEdgeData)>,
    }

    pub(super) fn arb_graph(alphabet: AnnotationAlphabet) -> impl Strategy<Value = NewickGraph> {
      let leaf = arb_node(alphabet, false).prop_map(|node| GenNode {
        node,
        children: Vec::new(),
      });
      let tree = leaf.prop_recursive(4, 40, 4, move |inner| {
        (arb_node(alphabet, true), vec((inner, arb_edge(alphabet)), 1..4))
          .prop_map(|(node, children)| GenNode { node, children })
      });
      tree.prop_map(|root| {
        let mut graph = NewickGraph::new();
        graph.root = add_subtree(&mut graph, root);
        graph
      })
    }

    pub(super) fn without_annotations(graph: &NewickGraph) -> NewickGraph {
      let mut graph = graph.clone();
      for node in &mut graph.nodes {
        node.node_attrs.clear();
        node.raw_comments.clear();
      }
      for edge in &mut graph.edges {
        edge.data.branch_attrs.clear();
        edge.data.raw_comments.clear();
      }
      graph
    }

    fn add_subtree(graph: &mut NewickGraph, node: GenNode) -> usize {
      let children: Vec<(usize, NewickEdgeData)> = node
        .children
        .into_iter()
        .map(|(child, edge)| (add_subtree(graph, child), edge))
        .collect();
      let idx = graph.add_node(node.node);
      for (child, edge) in children {
        graph.add_edge(idx, child, edge);
      }
      idx
    }

    fn arb_node(alphabet: AnnotationAlphabet, is_internal: bool) -> impl Strategy<Value = NewickNodeData> {
      (arb_label(is_internal), arb_attrs(alphabet), arb_raw_comments()).prop_map(|(label, node_attrs, raw_comments)| {
        NewickNodeData {
          label,
          node_attrs,
          raw_comments,
          hybrid: None,
          children: Vec::new(),
        }
      })
    }

    fn arb_edge(alphabet: AnnotationAlphabet) -> impl Strategy<Value = NewickEdgeData> {
      (
        proptest::option::of(arb_finite()),
        arb_attrs(alphabet),
        arb_raw_comments(),
      )
        .prop_map(|(branch_length, branch_attrs, raw_comments)| NewickEdgeData {
          branch_length,
          branch_attrs,
          raw_comments,
          is_acceptor: false,
        })
    }

    fn arb_label(is_internal: bool) -> BoxedStrategy<Option<NewickLabel>> {
      let any_name = "\\PC{0,8}";
      if is_internal {
        prop_oneof![
          Just(None),
          arb_finite().prop_map(|support| Some(NewickLabel::Support(support))),
          any_name
            .prop_filter("an internal name that parses as a number is a support value", |name| {
              name.parse::<f64>().is_err()
            })
            .prop_map(|name| Some(NewickLabel::Name(name))),
        ]
        .boxed()
      } else {
        prop_oneof![Just(None), any_name.prop_map(|name| Some(NewickLabel::Name(name)))].boxed()
      }
    }

    fn arb_attrs(alphabet: AnnotationAlphabet) -> BoxedStrategy<BTreeMap<String, NewickValue>> {
      match alphabet {
        AnnotationAlphabet::Beast => btree_map("\\PC{1,6}", arb_beast_value(), 0..3).boxed(),
        AnnotationAlphabet::Nhx => btree_map(
          "[A-Za-z0-9_]{1,6}",
          "[^:=\\[\\]\"\\s]{1,6}".prop_map(NewickValue::String),
          0..3,
        )
        .boxed(),
      }
    }

    fn arb_beast_value() -> impl Strategy<Value = NewickValue> {
      let scalar = prop_oneof![
        any::<bool>().prop_map(NewickValue::Boolean),
        arb_finite().prop_map(NewickValue::Number),
        "[0-9]{1,3}\\.[0-9]{0,2}0".prop_map(NewickValue::NumberText),
        "\\PC{0,8}".prop_map(NewickValue::String),
      ];
      scalar.prop_recursive(2, 8, 3, |inner| vec(inner, 0..3).prop_map(NewickValue::Array))
    }

    fn arb_raw_comments() -> impl Strategy<Value = Vec<String>> {
      vec("\\[[a-z ]{0,6}(\\[[a-z]{0,3}\\])?\\]", 0..2)
    }

    fn arb_finite() -> impl Strategy<Value = f64> {
      proptest::num::f64::NORMAL | proptest::num::f64::SUBNORMAL | proptest::num::f64::ZERO
    }
  }
}
