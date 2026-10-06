#[cfg(test)]
mod tests {
  use crate::__tests__::test_read_basic::tests::helpers::{caterpillar, on_small_stack, read_with};
  use crate::model::comment::{EdgeComment, EdgeField, LabelSide, NewickComment, NodeComment, ValueSide};
  use crate::model::data::{NewickEdgeData, NewickHybrid, NewickNodeData, SupportSource};
  use crate::model::graph::NewickGraph;
  use crate::model::value::NewickValue;
  use crate::read::options::NewickReadOptions;
  use helpers::{full_graph, network};
  use pretty_assertions::assert_eq;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::unknown_node(    (0, 9),  "Cannot connect node 0 to node 9: the graph has 3 nodes")]
  #[case::into_root(       (1, 0),  "The root, node 0 ('R'), cannot have a parent")]
  #[case::self_loop(       (1, 1),  "node 1 ('A') cannot be its own parent")]
  #[case::duplicate(       (0, 1),  "node 0 ('R') already has an edge to node 1 ('A')")]
  #[case::second_parent(   (2, 1),  "node 1 ('A') already has a parent, and only a hybrid node can have more than one")]
  #[trace]
  fn test_model_add_edge_rejects(#[case] (parent, child): (usize, usize), #[case] expected: &str) {
    let mut graph = NewickGraph::new(NewickNodeData::new().with_name("R"));
    graph.add_child(0, NewickEdgeData::new(), NewickNodeData::new().with_name("A")).unwrap();
    graph.add_child(0, NewickEdgeData::new(), NewickNodeData::new().with_name("B")).unwrap();

    let actual = graph.add_edge(parent, child, NewickEdgeData::new()).unwrap_err();

    assert_eq!(expected, actual.to_string());
  }

  #[test]
  fn test_model_hybrid_takes_second_parent() {
    let graph = network();

    assert_eq!(
      (vec![0, 2], vec![0, 3]),
      (graph.parents(1).collect::<Vec<_>>(), graph.parent_edges(1).to_vec())
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::preorder(   true,  vec![0, 1, 4, 2, 3])]
  #[case::postorder(  false, vec![4, 1, 3, 2, 0])]
  #[trace]
  fn test_model_traversal_visits_each_node_once(#[case] is_preorder: bool, #[case] expected: Vec<usize>) {
    let graph = network();

    let actual: Vec<usize> = if is_preorder { graph.preorder().collect() } else { graph.postorder().collect() };

    assert_eq!(expected, actual);
  }

  #[test]
  fn test_model_validate_rejects_unreachable_node() {
    let mut graph = NewickGraph::new(NewickNodeData::new());
    graph.add_node(NewickNodeData::new().with_name("lost"));

    assert_eq!(
      "node 1 ('lost') is not reachable from the root",
      graph.validate().unwrap_err().to_string()
    );
  }

  #[test]
  fn test_model_validate_rejects_cycle() {
    let mut graph = NewickGraph::new(NewickNodeData::new());
    let hybrid = graph
      .add_child(
        0,
        NewickEdgeData::new(),
        NewickNodeData::new().with_hybrid(NewickHybrid::new(None, 1)),
      )
      .unwrap();
    let below = graph
      .add_child(hybrid, NewickEdgeData::new(), NewickNodeData::new())
      .unwrap();
    graph.add_edge(below, hybrid, NewickEdgeData::new()).unwrap();

    assert_eq!("The graph contains a cycle", graph.validate().unwrap_err().to_string());
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::rooted(            |graph: &mut NewickGraph| graph.set_rooted(Some(false)))]
  #[case::weight(            |graph: &mut NewickGraph| graph.set_weight(Some(0.5)))]
  #[case::root_length(       |graph: &mut NewickGraph| graph.root_edge_mut().set_branch_length(Some(2.0)))]
  #[case::root_comment(      |graph: &mut NewickGraph| graph.root_edge_mut().comments_mut().clear())]
  #[case::name(              |graph: &mut NewickGraph| graph.node_mut(1).set_name(Some("Z".to_owned())))]
  #[case::hybrid(            |graph: &mut NewickGraph| *graph.node_mut(3) = NewickNodeData::new().with_name("C").with_hybrid(NewickHybrid::new(Some("LGT".to_owned()), 1)))]
  #[case::node_comment(      |graph: &mut NewickGraph| graph.node_mut(1).comments_mut()[0].comment = NewickComment::Plain("other".to_owned()))]
  #[case::comment_position(  |graph: &mut NewickGraph| graph.node_mut(1).comments_mut()[0].position = LabelSide::BeforeLabel)]
  #[case::branch_length(     |graph: &mut NewickGraph| graph.edge_mut(0).set_branch_length(Some(9.0)))]
  #[case::support(           |graph: &mut NewickGraph| graph.edge_mut(0).set_support(vec![81.0], SupportSource::Field))]
  #[case::support_source(    |graph: &mut NewickGraph| graph.edge_mut(0).set_support(vec![80.0], SupportSource::Label))]
  #[case::probability(       |graph: &mut NewickGraph| graph.edge_mut(0).set_probability(None))]
  #[case::edge_comment_side( |graph: &mut NewickGraph| graph.edge_mut(0).comments_mut()[0].side = ValueSide::AfterValue)]
  #[case::edge_comment_field(|graph: &mut NewickGraph| graph.edge_mut(0).comments_mut()[0].field = EdgeField::Support)]
  #[case::annotation_value(  |graph: &mut NewickGraph| graph.edge_mut(0).comments_mut()[0].comment = NewickComment::Beast(vec![("rate".to_owned(), NewickValue::NumberText("1.50".to_owned()))]))]
  #[case::acceptor(          |graph: &mut NewickGraph| *graph.edge_mut(2) = NewickEdgeData::new().with_acceptor(false))]
  fn test_model_equality_sees_every_field(#[case] change: fn(&mut NewickGraph)) {
    let original = full_graph();
    let mut changed = original.clone();

    change(&mut changed);

    assert_eq!((true, false), (original == original.clone(), original == changed));
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::reordered(         "((A,B),(C,D));",   "((D,C),(B,A));",   (true,  false))]
  #[case::other_topology(    "((A,B),C);",       "(A,(B,C));",       (false, false))]
  #[case::same_multiset(     "((A,A),(A,B));",   "((A,B),(A,A));",   (true,  false))]
  #[case::other_multiset(    "((A,A),(B,B));",   "((A,B),(A,B));",   (false, false))]
  #[case::same_order(        "(A:1,B:2);",       "(A:1,B:2);",       (true,  true))]
  #[trace]
  fn test_model_equality_unordered_and_ordered(#[case] left: &str, #[case] right: &str, #[case] expected: (bool, bool)) {
    let options = NewickReadOptions::default();
    let left = read_with(left, &options).graph;
    let right = read_with(right, &options).graph;

    assert_eq!(expected, (left == right, left.eq_ordered(&right)));
  }

  #[test]
  fn test_model_equality_deep_tree_on_small_stack() {
    let input = caterpillar(100_000);

    let equal = on_small_stack(move || {
      let graph = read_with(&input, &NewickReadOptions::default()).graph;
      graph == graph.clone()
    });

    assert!(equal);
  }

  #[test]
  fn test_model_annotations_of_edge() {
    let edge = NewickEdgeData::new()
      .with_comment(EdgeComment::new(
        EdgeField::Length,
        ValueSide::BeforeValue,
        NewickComment::Beast(vec![("a".to_owned(), NewickValue::Number(1.0))]),
      ))
      .with_comment(EdgeComment::new(
        EdgeField::Length,
        ValueSide::AfterValue,
        NewickComment::Nhx(vec![("a".to_owned(), NewickValue::Number(2.0))]),
      ));

    assert_eq!(
      (2, Some(&NewickValue::Number(2.0))),
      (edge.annotations().count(), edge.annotation("a"))
    );
  }

  #[test]
  fn test_model_node_comment_builder() {
    let node = NewickNodeData::new().with_comment(NodeComment::new(
      LabelSide::AfterLabel,
      NewickComment::Plain("x".to_owned()),
    ));

    assert_eq!(
      vec![NodeComment::new(
        LabelSide::AfterLabel,
        NewickComment::Plain("x".to_owned())
      )],
      node.comments().to_vec()
    );
  }

  mod helpers {
    use crate::model::comment::{EdgeComment, EdgeField, LabelSide, NewickComment, NodeComment, ValueSide};
    use crate::model::data::{NewickEdgeData, NewickHybrid, NewickNodeData, SupportSource};
    use crate::model::graph::NewickGraph;
    use crate::model::value::NewickValue;

    pub(super) fn network() -> NewickGraph {
      let hybrid = NewickNodeData::new().with_hybrid(NewickHybrid::new(Some("H".to_owned()), 1));
      let mut graph = NewickGraph::new(NewickNodeData::new());
      let h = graph.add_child(0, NewickEdgeData::new(), hybrid).unwrap();
      let inner = graph
        .add_child(0, NewickEdgeData::new(), NewickNodeData::new())
        .unwrap();
      graph
        .add_child(inner, NewickEdgeData::new(), NewickNodeData::new().with_name("B"))
        .unwrap();
      graph.add_edge(inner, h, NewickEdgeData::new()).unwrap();
      graph
        .add_child(h, NewickEdgeData::new(), NewickNodeData::new().with_name("A"))
        .unwrap();
      graph
    }

    pub(super) fn full_graph() -> NewickGraph {
      let mut graph = NewickGraph::new(NewickNodeData::new().with_name("R"));
      graph.set_rooted(Some(true));
      graph.set_weight(Some(1.0));
      *graph.root_edge_mut() = NewickEdgeData::new().with_length(1.0).with_comment(EdgeComment::new(
        EdgeField::Length,
        ValueSide::AfterValue,
        NewickComment::Plain("r".to_owned()),
      ));
      let leaf = NewickNodeData::new().with_name("A").with_comment(NodeComment::new(
        LabelSide::AfterLabel,
        NewickComment::Plain("a".to_owned()),
      ));
      let edge = NewickEdgeData::new()
        .with_length(1.0)
        .with_support(vec![80.0], SupportSource::Field)
        .with_probability(0.4)
        .with_comment(EdgeComment::new(
          EdgeField::Length,
          ValueSide::BeforeValue,
          NewickComment::Beast(vec![("rate".to_owned(), NewickValue::Number(1.5))]),
        ));
      graph.add_child(0, edge, leaf).unwrap();
      let inner = graph
        .add_child(0, NewickEdgeData::new(), NewickNodeData::new())
        .unwrap();
      let hybrid = NewickNodeData::new()
        .with_name("C")
        .with_hybrid(NewickHybrid::new(Some("H".to_owned()), 1));
      let h = graph
        .add_child(inner, NewickEdgeData::new().with_acceptor(true), hybrid)
        .unwrap();
      graph.add_edge(0, h, NewickEdgeData::new()).unwrap();
      graph
        .add_child(h, NewickEdgeData::new(), NewickNodeData::new().with_name("D"))
        .unwrap();
      graph
    }
  }
}
