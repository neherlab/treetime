#[cfg(test)]
mod tests {
  use super::super::test_graph_support::tests::helpers::NamedGraph;
  use crate::schedule::GraphSchedule;
  use crate::scope::{EdgeScopes, SiteScope};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use smallvec::smallvec;
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
}
