#[cfg(test)]
mod tests {
  use self::helpers::TopologyChange;
  use super::super::test_graph_support::tests::helpers::NamedGraph;
  use crate::schedule::{GraphSchedule, ScheduleIndex};
  use eyre::Report;
  use pretty_assertions::{assert_eq, assert_ne};
  use rstest::rstest;

  #[rustfmt::skip]
  #[rstest]
  #[case::add_leaf(                 TopologyChange::AddLeaf)]
  #[case::add_root(                 TopologyChange::AddRoot)]
  #[case::remove_edge(              TopologyChange::RemoveEdge)]
  #[case::collapse_edge(            TopologyChange::CollapseEdge)]
  #[case::reparent_edge(            TopologyChange::ReparentEdge)]
  #[case::invert_edge(              TopologyChange::InvertEdge)]
  #[case::order_by_descendant_count(TopologyChange::OrderByDescendantCount)]
  #[case::rebuild_roots(            TopologyChange::RebuildRoots)]
  #[trace]
  fn test_schedule_cache_follows_each_topology_change(#[case] change: TopologyChange) -> Result<(), Report> {
    let mut tree = NamedGraph::new(
      &["root", "AB", "A", "B", "C", "X"],
      &[("root", "AB"), ("root", "C"), ("AB", "A"), ("AB", "B")],
    )?;
    let before = ScheduleIndex::build(&tree.graph)?;
    GraphSchedule::new(&tree.graph)?;

    change.apply(&mut tree)?;

    let expected = ScheduleIndex::build(&tree.graph)?;
    assert_ne!(before, expected);
    assert_eq!(&expected, tree.graph.schedule_index()?);
    Ok(())
  }

  mod helpers {
    use super::super::super::test_graph_support::tests::helpers::NamedGraph;
    use crate::edge::invert_edge;
    use crate::node_names::NodeNames;
    use crate::schedule::GraphSchedule;
    use crate::topology_order::TopologyOrderSpec;
    use eyre::Report;
    use std::collections::BTreeMap;

    #[derive(Clone, Copy, Debug)]
    pub(super) enum TopologyChange {
      AddLeaf,
      AddRoot,
      RemoveEdge,
      CollapseEdge,
      ReparentEdge,
      InvertEdge,
      OrderByDescendantCount,
      RebuildRoots,
    }

    impl TopologyChange {
      pub(super) fn apply(self, tree: &mut NamedGraph) -> Result<(), Report> {
        match self {
          Self::AddLeaf => {
            let leaf = tree.graph.add_node();
            tree.graph.add_edge(tree.key("C"), leaf)?;
          },
          Self::AddRoot => {
            tree.graph.add_node();
            tree.graph.build()?;
          },
          Self::RemoveEdge => {
            let edge = tree.edge("root", "C")?;
            tree.graph.remove_edge(edge)?;
          },
          Self::CollapseEdge => {
            let edge = tree.edge("root", "AB")?;
            tree.graph.collapse_edge(edge)?;
          },
          Self::ReparentEdge => {
            let (edge, source) = (tree.edge("AB", "B")?, tree.key("C"));
            tree.graph.reparent_edge(edge, source)?;
          },
          Self::InvertEdge => {
            let edge = tree.edge("root", "AB")?;
            invert_edge(&mut tree.graph, edge);
          },
          Self::OrderByDescendantCount => order_by_descendant_count(tree)?,
          Self::RebuildRoots => {
            order_by_descendant_count(tree)?;
            GraphSchedule::new(&tree.graph)?;
            tree.graph.build()?;
          },
        }
        Ok(())
      }
    }

    fn order_by_descendant_count(tree: &mut NamedGraph) -> Result<(), Report> {
      let names = NodeNames::from_input(BTreeMap::new(), &tree.graph)?.names;
      TopologyOrderSpec::default().apply(&mut tree.graph, &names, &BTreeMap::new())
    }
  }
}
