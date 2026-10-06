#[cfg(test)]
mod tests {
  use crate::dialect::NewickDialect;
  use crate::number::NumberFormat;
  use crate::read::options::NewickReadOptions;
  use crate::read::stream::newick_from_str;
  use crate::write::newick::newick_to_string;
  use crate::write::options::{NewickWriteOptions, Quoting};
  use generators::gen_graph;
  use proptest::prelude::*;
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::classic(  NewickDialect::Classic)]
  #[case::beast(    NewickDialect::Beast)]
  #[case::mrbayes(  NewickDialect::MrBayes)]
  #[case::nhx(      NewickDialect::Nhx)]
  #[case::enewick(  NewickDialect::ENewick)]
  #[case::rich(     NewickDialect::Rich)]
  #[trace]
  fn test_prop_roundtrip_dialect(#[case] dialect: NewickDialect) {
    proptest!(|(graph in gen_graph(dialect), quoting in prop_oneof![Just(Quoting::WhenNeeded), Just(Quoting::Always)], indent in proptest::option::of(0_usize..3), point_zero in any::<bool>())| {
      let options = NewickWriteOptions {
        quoting,
        indent,
        numbers: NumberFormat { point_zero, ..NumberFormat::default() },
        ..NewickWriteOptions::new(dialect)
      };
      let read_options = NewickReadOptions { dialects: vec![dialect], ..NewickReadOptions::default() };

      let written = newick_to_string(&graph, &options).unwrap();
      let tree = newick_from_str(&written, &read_options).map_err(|error| TestCaseError::fail(format!("{error}\nWritten: {written}")))?;

      prop_assert!(graph.eq_ordered(&tree.graph), "Round trip changed the graph.\nWritten: {written}\nBefore: {graph:#?}\nAfter: {:#?}", tree.graph);
    });
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::beast(    NewickDialect::Beast)]
  #[case::rich(     NewickDialect::Rich)]
  #[trace]
  fn test_prop_roundtrip_write_idempotent(#[case] dialect: NewickDialect) {
    proptest!(|(graph in gen_graph(dialect))| {
      let options = NewickWriteOptions::new(dialect);
      let read_options = NewickReadOptions { dialects: vec![dialect], ..NewickReadOptions::default() };

      let first = newick_to_string(&graph, &options).unwrap();
      let second = newick_to_string(&newick_from_str(&first, &read_options).unwrap().graph, &options).unwrap();

      prop_assert_eq!(first, second);
    });
  }

  #[test]
  fn test_prop_roundtrip_all_dialects_read_back_rich_output() {
    proptest!(|(graph in gen_graph(NewickDialect::Rich))| {
      let written = newick_to_string(&graph, &NewickWriteOptions::new(NewickDialect::Rich)).unwrap();

      let tree = newick_from_str(&written, &NewickReadOptions::all_dialects()).unwrap();

      prop_assert!(graph.eq_ordered(&tree.graph), "Written: {written}");
    });
  }

  mod generators {
    use crate::dialect::{CommentKind, NewickDialect};
    use crate::model::comment::{
      EdgeComment, EdgeField, LabelSide, MrBayesComment, MrBayesKind, NewickComment, NodeComment, ValueSide,
    };
    use crate::model::data::{NewickEdgeData, NewickHybrid, NewickNodeData, SupportSource};
    use crate::model::graph::NewickGraph;
    use crate::model::value::NewickValue;
    use proptest::collection::vec;
    use proptest::prelude::*;

    #[derive(Clone, Debug)]
    pub(super) struct GenNode {
      name: Option<String>,
      hybrid: bool,
      comments: Vec<(bool, NewickComment)>,
      children: Vec<(GenNode, GenEdge)>,
    }

    #[derive(Clone, Debug)]
    pub(super) struct GenEdge {
      length: Option<f64>,
      support: Option<(Vec<f64>, bool)>,
      probability: Option<f64>,
      acceptor: bool,
      comments: Vec<(usize, bool, NewickComment)>,
    }

    pub(super) fn gen_graph(dialect: NewickDialect) -> impl Strategy<Value = NewickGraph> {
      let leaf = gen_node_fields(dialect).prop_map(|(name, hybrid, comments)| GenNode {
        name,
        hybrid,
        comments,
        children: Vec::new(),
      });
      let tree = leaf.prop_recursive(4, 32, 4, move |inner| {
        (gen_node_fields(dialect), vec((inner, gen_edge(dialect)), 1..4)).prop_map(
          |((name, hybrid, comments), children)| GenNode {
            name,
            hybrid,
            comments,
            children,
          },
        )
      });
      let extra = vec((any::<usize>(), any::<usize>(), gen_edge(dialect)), 0..3);
      let tree_comments = (proptest::option::of(any::<bool>()), proptest::option::of(gen_finite()));
      (tree, gen_edge(dialect), extra, tree_comments).prop_map(move |(root, root_edge, extra, (rooted, weight))| {
        Builder::new(dialect).build(&root, &root_edge, &extra, rooted, weight)
      })
    }

    struct Builder {
      dialect: NewickDialect,
      next_hybrid: u32,
    }

    impl Builder {
      fn new(dialect: NewickDialect) -> Self {
        Self {
          dialect,
          next_hybrid: 0,
        }
      }

      fn build(
        mut self,
        root: &GenNode,
        root_edge: &GenEdge,
        extra: &[(usize, usize, GenEdge)],
        rooted: Option<bool>,
        weight: Option<f64>,
      ) -> NewickGraph {
        let features = self.dialect.features();
        let mut graph = NewickGraph::new(self.node_data(root, false));
        let is_internal = !root.children.is_empty();
        *graph.root_edge_mut() = self.edge_data(root_edge, is_internal, root.name.is_some(), false);
        self.add_children(&mut graph, 0, root);
        let hybrids: Vec<usize> = graph
          .nodes()
          .filter(|(_, node)| node.hybrid().is_some())
          .map(|(idx, _)| idx)
          .collect();
        for (hybrid_pick, parent_pick, edge) in extra {
          let Some(&hybrid) = hybrids.get(hybrid_pick % hybrids.len().max(1)) else {
            continue;
          };
          let below = reachable_from(&graph, hybrid);
          let candidates: Vec<usize> = (0..graph.node_count())
            .filter(|&node| !below[node] && !graph.children(node).any(|child| child == hybrid))
            .collect();
          if let Some(&parent) = candidates.get(parent_pick % candidates.len().max(1)) {
            let data = self.edge_data(edge, false, true, true);
            graph.add_edge(parent, hybrid, data).unwrap();
          }
        }
        graph.set_rooted(rooted.filter(|_| features.rooting));
        graph.set_weight(weight.filter(|_| features.weight));
        graph
      }

      fn add_children(&mut self, graph: &mut NewickGraph, parent: usize, node: &GenNode) {
        for (child, edge) in &node.children {
          let data = self.node_data(child, true);
          let has_label = data.name().is_some() || data.hybrid().is_some();
          let is_hybrid = data.hybrid().is_some();
          let edge = self.edge_data(edge, !child.children.is_empty(), has_label, is_hybrid);
          let idx = graph.add_child(parent, edge, data).unwrap();
          self.add_children(graph, idx, child);
        }
      }

      fn node_data(&mut self, node: &GenNode, may_be_hybrid: bool) -> NewickNodeData {
        let mut data = NewickNodeData::new();
        if let Some(name) = &node.name {
          data = data.with_name(name.clone());
        }
        if node.hybrid && may_be_hybrid && self.dialect.features().hybrid_tags {
          self.next_hybrid += 1;
          data = data.with_hybrid(NewickHybrid::new(Some("H".to_owned()), self.next_hybrid));
        }
        let has_label = data.name().is_some() || data.hybrid().is_some();
        for (before, comment) in &node.comments {
          let position = if *before && has_label {
            LabelSide::BeforeLabel
          } else {
            LabelSide::AfterLabel
          };
          data = data.with_comment(NodeComment::new(position, comment.clone()));
        }
        data
      }

      fn edge_data(
        &self,
        edge: &GenEdge,
        child_internal: bool,
        child_labeled: bool,
        child_hybrid: bool,
      ) -> NewickEdgeData {
        let features = self.dialect.features();
        let mut data = NewickEdgeData::new().with_acceptor(edge.acceptor && child_hybrid);
        data.set_branch_length(edge.length);
        match &edge.support {
          Some((values, true)) if features.rich_fields && !values.is_empty() => {
            data.set_support(vec![values[0]], SupportSource::Field);
          },
          Some((values, false)) if child_internal && !child_labeled && !child_hybrid && !values.is_empty() => {
            data.set_support(values.clone(), SupportSource::Label);
          },
          Some(_) | None => {},
        }
        if features.rich_fields {
          data.set_probability(edge.probability);
        }
        let has_value = [
          data.branch_length().is_some(),
          data.support_source() == SupportSource::Field && !data.support().is_empty(),
          data.probability().is_some(),
        ];
        for (field, after, comment) in &edge.comments {
          let index = field % features.field_count();
          let side = if *after && has_value[index] {
            ValueSide::AfterValue
          } else {
            ValueSide::BeforeValue
          };
          data = data.with_comment(EdgeComment::new(EdgeField::ALL[index], side, comment.clone()));
        }
        data
      }
    }

    fn reachable_from(graph: &NewickGraph, start: usize) -> Vec<bool> {
      let mut reached = vec![false; graph.node_count()];
      let mut pending = vec![start];
      while let Some(node) = pending.pop() {
        if !std::mem::replace(&mut reached[node], true) {
          pending.extend(graph.children(node));
        }
      }
      reached
    }

    fn gen_node_fields(
      dialect: NewickDialect,
    ) -> impl Strategy<Value = (Option<String>, bool, Vec<(bool, NewickComment)>)> {
      (
        proptest::option::of("\\PC{0,8}"),
        proptest::bool::weighted(0.2),
        vec((any::<bool>(), gen_comment(dialect)), 0..2),
      )
    }

    fn gen_edge(dialect: NewickDialect) -> impl Strategy<Value = GenEdge> {
      (
        proptest::option::of(gen_finite()),
        proptest::option::of((vec(gen_finite(), 1..3), any::<bool>())),
        proptest::option::of(gen_finite()),
        any::<bool>(),
        vec((0_usize..3, any::<bool>(), gen_comment(dialect)), 0..2),
      )
        .prop_map(|(length, support, probability, acceptor, comments)| GenEdge {
          length,
          support,
          probability,
          acceptor,
          comments,
        })
    }

    fn gen_comment(dialect: NewickDialect) -> BoxedStrategy<NewickComment> {
      let features = dialect.features();
      let plain = "[a-z \"'=,{}:;()#&]{0,6}(\\[[a-z;]{0,3}\\])?[a-z ]{0,3}"
        .prop_filter(
          "a reserved dialect reads '[&' as an annotation",
          move |text: &String| !(features.reserves_annotations && text.starts_with('&')),
        )
        .prop_map(NewickComment::Plain);
      match features.comments {
        CommentKind::Plain => plain.boxed(),
        CommentKind::Beast => prop_oneof![plain, gen_beast_comment()].boxed(),
        CommentKind::Nhx => prop_oneof![plain, gen_nhx_comment()].boxed(),
        CommentKind::MrBayes => prop_oneof![plain, gen_mrbayes_comment()].boxed(),
      }
    }

    fn gen_beast_comment() -> impl Strategy<Value = NewickComment> {
      vec(("\\PC{0,6}", gen_beast_value()), 0..3).prop_map(NewickComment::Beast)
    }

    fn gen_beast_value() -> impl Strategy<Value = NewickValue> {
      let scalar = prop_oneof![
        any::<bool>().prop_map(NewickValue::Boolean),
        gen_finite().prop_map(NewickValue::Number),
        "[0-9]{1,3}\\.[0-9]{0,2}0".prop_map(NewickValue::NumberText),
        "\\PC{0,8}".prop_map(NewickValue::String),
        any::<[u8; 3]>().prop_map(NewickValue::Color),
      ];
      scalar.prop_recursive(2, 8, 3, |inner| {
        vec(inner, 0..3).prop_map(|values| NewickValue::Array(values.into()))
      })
    }

    fn gen_nhx_comment() -> impl Strategy<Value = NewickComment> {
      let text_part = "[^:=\\[\\]>]{0,6}";
      let tag = prop_oneof![
        ("[a-z][a-z0-9_]{0,4}", text_part).prop_map(|(key, value)| (key, NewickValue::String(value))),
        ("[a-z][a-z0-9_]{0,4}", vec(text_part, 2..4)).prop_map(|(key, parts)| {
          (
            key,
            NewickValue::Array(parts.into_iter().map(NewickValue::String).collect::<Vec<_>>().into()),
          )
        }),
        "[a-z][a-z0-9_]{0,4}".prop_map(|key| (key, NewickValue::Boolean(true))),
        gen_finite().prop_map(|value| ("B".to_owned(), NewickValue::Number(value))),
        (-1000_i32..1000).prop_map(|value| ("T".to_owned(), NewickValue::Number(f64::from(value)))),
        prop_oneof![Just("T"), Just("F"), Just("Y"), Just("N"), Just("?")]
          .prop_map(|value| ("D".to_owned(), NewickValue::String(value.to_owned()))),
        any::<[u8; 3]>().prop_map(|value| ("C".to_owned(), NewickValue::Color(value))),
      ];
      vec(tag, 0..3).prop_map(NewickComment::Nhx)
    }

    fn gen_mrbayes_comment() -> impl Strategy<Value = NewickComment> {
      let token = "[A-Za-z0-9_.:=]{1,6}";
      (
        prop_oneof![Just(MrBayesKind::E), Just(MrBayesKind::B), Just(MrBayesKind::N)],
        token,
        vec(token, 0..3),
      )
        .prop_map(|(kind, name, values)| NewickComment::MrBayesMcmc(MrBayesComment { kind, name, values }))
    }

    fn gen_finite() -> impl Strategy<Value = f64> {
      proptest::num::f64::NORMAL | proptest::num::f64::SUBNORMAL | proptest::num::f64::ZERO
    }
  }
}
