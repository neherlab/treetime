#[cfg(test)]
mod tests {
  use super::super::test_graph_support::tests::helpers::NamedGraph;
  use crate::schedule::GraphSchedule;
  use crate::scope::{EdgeScopes, Inheritance, SiteScope};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use smallvec::smallvec;
  use std::sync::Arc;
  use treetime_utils::assert_error;

  fn ranges(ranges: &[(usize, usize)]) -> SiteScope {
    SiteScope::from_ranges(ranges.iter().copied())
  }

  #[test]
  fn test_site_scope_from_ranges_sorts_merges_and_drops_empty_ranges() {
    assert_eq!(
      SiteScope::Ranges(smallvec![(0, 10), (20, 30)]),
      ranges(&[(25, 30), (5, 10), (7, 7), (0, 6), (20, 25)])
    );
  }

  #[test]
  fn test_site_scope_from_no_ranges_is_empty() {
    assert_eq!(SiteScope::Empty, ranges(&[(3, 3)]));
  }

  #[rstest]
  #[case(SiteScope::Empty, 0, false)]
  #[case(SiteScope::Full, 1_000_000, true)]
  #[case(ranges(&[(10, 20)]), 9, false)]
  #[case(ranges(&[(10, 20)]), 10, true)]
  #[case(ranges(&[(10, 20)]), 19, true)]
  #[case(ranges(&[(10, 20)]), 20, false)]
  fn test_site_scope_contains_half_open_ranges(#[case] scope: SiteScope, #[case] site: usize, #[case] expected: bool) {
    assert_eq!(expected, scope.contains(site));
  }

  #[rstest]
  #[case(SiteScope::Full, SiteScope::Empty, false)]
  #[case(SiteScope::Full, ranges(&[(0, 1)]), true)]
  #[case(ranges(&[(0, 10)]), ranges(&[(10, 20)]), false)]
  #[case(ranges(&[(0, 11)]), ranges(&[(10, 20)]), true)]
  fn test_site_scope_overlaps(#[case] a: SiteScope, #[case] b: SiteScope, #[case] expected: bool) {
    assert_eq!((expected, expected), (a.overlaps(&b), b.overlaps(&a)));
  }

  #[rstest]
  #[case(SiteScope::Empty, SiteScope::Empty, true)]
  #[case(SiteScope::Full, SiteScope::Full, true)]
  #[case(SiteScope::Full, ranges(&[(0, 100)]), false)]
  #[case(ranges(&[(0, 100)]), SiteScope::Full, true)]
  #[case(ranges(&[(5, 10)]), ranges(&[(0, 10), (20, 30)]), true)]
  #[case(ranges(&[(5, 21)]), ranges(&[(0, 10), (20, 30)]), false)]
  #[case(ranges(&[(5, 10)]), SiteScope::Empty, false)]
  fn test_site_scope_is_subset_of(#[case] a: SiteScope, #[case] b: SiteScope, #[case] expected: bool) {
    assert_eq!(expected, a.is_subset_of(&b));
  }

  #[test]
  fn test_site_scope_union_of_crossover_halves_is_one_range() {
    assert_eq!(ranges(&[(0, 100)]), ranges(&[(0, 40)]).union(&ranges(&[(40, 100)])));
  }

  #[test]
  fn test_edge_scopes_edge_without_entry_carries_nothing() -> Result<(), Report> {
    let tree = NamedGraph::new(&["r", "a"], &[("r", "a")])?;
    let scopes = EdgeScopes::new(2);
    assert_eq!(
      (false, &SiteScope::Empty),
      (
        scopes.carries(tree.edge("r", "a")?, 1),
        scopes.get(tree.edge("r", "a")?, 1)
      )
    );
    Ok(())
  }

  #[test]
  fn test_edge_scopes_set_rejects_segment_out_of_range() -> Result<(), Report> {
    let tree = NamedGraph::new(&["r", "a"], &[("r", "a")])?;
    let mut scopes = EdgeScopes::new(2);
    assert_error!(
      scopes.set(tree.edge("r", "a")?, 2, SiteScope::Full),
      "Segment 2 is out of range: the edge scopes have 2 segments. This is an internal error. Please report it to developers."
    );
    Ok(())
  }

  fn recombinant() -> Result<NamedGraph, Report> {
    NamedGraph::new(
      &["r", "p", "q", "c", "x", "y"],
      &[("r", "p"), ("r", "q"), ("p", "c"), ("q", "c"), ("c", "x"), ("c", "y")],
    )
  }

  fn scopes_of(graph: &NamedGraph, edges: &[(&str, &str, SiteScope)]) -> Result<EdgeScopes, Report> {
    let mut scopes = EdgeScopes::new(1);
    for (source, target, scope) in edges {
      scopes.set(graph.edge(source, target)?, 0, scope.clone())?;
    }
    Ok(scopes)
  }

  fn validate(graph: &NamedGraph, scopes: &EdgeScopes) -> Result<(), Report> {
    scopes.validate(&GraphSchedule::new(&graph.graph)?, |key| graph.name(key).to_owned())
  }

  #[test]
  fn test_edge_scopes_validate_accepts_crossover_with_disjoint_parent_ranges() -> Result<(), Report> {
    let graph = recombinant()?;
    let scopes = scopes_of(
      &graph,
      &[
        ("r", "p", ranges(&[(0, 40)])),
        ("r", "q", ranges(&[(40, 100)])),
        ("p", "c", ranges(&[(0, 40)])),
        ("q", "c", ranges(&[(40, 100)])),
        ("c", "x", SiteScope::Full),
        ("c", "y", SiteScope::Full),
      ],
    )?;
    validate(&graph, &scopes)
  }

  #[test]
  fn test_edge_scopes_validate_rejects_a_site_with_two_parents() -> Result<(), Report> {
    let graph = recombinant()?;
    let scopes = scopes_of(
      &graph,
      &[
        ("r", "p", ranges(&[(0, 41)])),
        ("r", "q", ranges(&[(40, 100)])),
        ("p", "c", ranges(&[(0, 41)])),
        ("q", "c", ranges(&[(40, 100)])),
        ("c", "x", SiteScope::Full),
        ("c", "y", SiteScope::Full),
      ],
    )?;
    assert_error!(
      validate(&graph, &scopes),
      "Node 'c' inherits sites of segment 1 from two parents, 'p' and 'q', but every site has one parent"
    );
    Ok(())
  }

  #[test]
  fn test_edge_scopes_validate_rejects_sites_that_no_child_edge_carries() -> Result<(), Report> {
    let graph = recombinant()?;
    let scopes = scopes_of(
      &graph,
      &[
        ("r", "p", ranges(&[(0, 50)])),
        ("r", "q", ranges(&[(40, 100)])),
        ("p", "c", ranges(&[(0, 40)])),
        ("q", "c", ranges(&[(40, 100)])),
        ("c", "x", SiteScope::Full),
        ("c", "y", SiteScope::Full),
      ],
    )?;
    assert_error!(
      validate(&graph, &scopes),
      "The edge 'r' -> 'p' carries sites of segment 1 that no child edge of 'p' carries"
    );
    Ok(())
  }

  #[test]
  fn test_edge_scopes_validate_accepts_a_full_tree() -> Result<(), Report> {
    let graph = NamedGraph::new(&["r", "a", "b"], &[("r", "a"), ("r", "b")])?;
    let scopes = scopes_of(&graph, &[("r", "a", SiteScope::Full), ("r", "b", SiteScope::Full)])?;
    validate(&graph, &scopes)
  }

  fn reassortment() -> Result<(NamedGraph, Arc<EdgeScopes>), Report> {
    let graph = NamedGraph::new(
      &["r1", "r2", "p", "c", "x", "y"],
      &[("r1", "p"), ("r2", "c"), ("p", "c"), ("c", "x"), ("c", "y")],
    )?;
    let mut scopes = EdgeScopes::new(2);
    for (source, target, segments) in [
      ("r1", "p", vec![0]),
      ("p", "c", vec![0]),
      ("r2", "c", vec![1]),
      ("c", "x", vec![0, 1]),
      ("c", "y", vec![0, 1]),
    ] {
      for segment in segments {
        scopes.set(graph.edge(source, target)?, segment, SiteScope::Full)?;
      }
    }
    Ok((graph, Arc::new(scopes)))
  }

  #[test]
  fn test_inheritance_follows_the_parent_that_carries_each_segment() -> Result<(), Report> {
    let (graph, scopes) = reassortment()?;
    let schedule = GraphSchedule::new(&graph.graph)?;
    let parent_of_c = |segment| -> Result<Option<&str>, Report> {
      let inheritance = Inheritance::segment(Arc::clone(&scopes), segment)?;
      Ok(
        inheritance
          .parent(&schedule, graph.key("c"))
          .map(|(_, parent, _)| graph.name(parent)),
      )
    };
    assert_eq!((Some("p"), Some("r2")), (parent_of_c(0)?, parent_of_c(1)?));
    Ok(())
  }

  #[test]
  fn test_inheritance_holds_only_the_nodes_on_edges_that_carry_the_segment() -> Result<(), Report> {
    let (graph, scopes) = reassortment()?;
    let schedule = GraphSchedule::new(&graph.graph)?;
    let inheritance = Inheritance::segment(scopes, 1)?;
    let held = schedule
      .preorder()
      .iter()
      .filter(|&&key| inheritance.holds(&schedule, key))
      .map(|&key| graph.name(key))
      .collect::<Vec<_>>();
    let children = inheritance
      .children(&schedule, graph.key("r2"))
      .map(|(_, child, _)| graph.name(child))
      .collect::<Vec<_>>();
    assert_eq!((vec!["r2", "c", "x", "y"], vec!["c"]), (held, children));
    Ok(())
  }

  #[test]
  fn test_inheritance_of_one_segment_accepts_a_node_with_two_parents_and_two_roots() -> Result<(), Report> {
    let (graph, scopes) = reassortment()?;
    let schedule = GraphSchedule::new(&graph.graph)?;
    for segment in 0..2 {
      Inheritance::segment(Arc::clone(&scopes), segment)?.require_one_parent(&schedule)?;
    }
    Ok(())
  }

  #[test]
  fn test_inheritance_of_all_sites_rejects_a_node_with_two_parents() -> Result<(), Report> {
    let (graph, _) = reassortment()?;
    let schedule = GraphSchedule::new(&graph.graph)?;
    let c = graph.key("c");
    assert_error!(
      Inheritance::ALL.require_one_parent(&schedule),
      format!("Node {c} has 2 parents that pass on all sites, but every site has one parent")
    );
    Ok(())
  }

  #[test]
  fn test_inheritance_rejects_an_edge_that_carries_part_of_a_segment() -> Result<(), Report> {
    let graph = NamedGraph::new(&["r", "a"], &[("r", "a")])?;
    let edge = graph.edge("r", "a")?;
    let mut scopes = EdgeScopes::new(1);
    scopes.set(edge, 0, ranges(&[(0, 10)]))?;
    let schedule = GraphSchedule::new(&graph.graph)?;
    assert_error!(
      Inheritance::segment(Arc::new(scopes), 0)?.require_one_parent(&schedule),
      format!(
        "Edge {edge} passes on part of segment 1, but this reconstruction supports only edges that carry a segment entirely or not at all"
      )
    );
    Ok(())
  }

  #[test]
  fn test_inheritance_of_a_segment_without_scopes_is_error() {
    assert_error!(
      Inheritance::segment(Arc::new(EdgeScopes::new(2)), 2),
      "Segment 3 has no edge scopes: the graph has scopes for 2 segments"
    );
  }
}
