#[cfg(test)]
mod tests {
  use crate::reroot::find_root::EdgePoint;
  use crate::reroot::placement::{Placement, RootTarget, leaf_keys, require_dated_new_leaves, root_moves};
  use crate::test_utils::find_node_key_by_name;
  use eyre::Report;
  use helpers::{Expected, edge_fixture, terminal_edge_fixture, zero_length_edge_fixture};
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use treetime_graph::split_fraction::SplitFraction;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::assert_error;

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(                         0.0,            true,  Expected::Source)]
  #[case::inside_parent_end_band(       9.9e-6,         true,  Expected::Source)]
  #[case::parent_end_band_boundary(     1e-5,           true,  Expected::Split)]
  #[case::one(                          1.0,            true,  Expected::Target)]
  #[case::inside_child_end_band(        1.0 - 9.9e-6,   true,  Expected::Target)]
  #[case::outside_child_end_band(       1.0 - 1.1e-5,   true,  Expected::Split)]
  #[case::interior_split(               0.25,           true,  Expected::Split)]
  #[case::interior_split_upper(         0.75,           true,  Expected::Split)]
  #[case::below_half_no_split(          0.25,           false, Expected::Source)]
  #[case::half_no_split(                0.5,            false, Expected::Target)]
  #[case::above_half_no_split(          0.75,           false, Expected::Target)]
  #[trace]
  fn test_root_target_on_edge(
    #[case] split: f64,
    #[case] split_edge: bool,
    #[case] expected: Expected,
  ) -> Result<(), Report> {
    let fixture = edge_fixture()?;
    let split = SplitFraction::new(split)?;
    let expected = expected.root_target(&fixture, split);

    let actual = fixture.target(split, Placement::EdgeEndBand, split_edge)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::zero(                    0.0,            true,  Expected::Source)]
  #[case::inside_parent_end_band(  9.9e-6,         true,  Expected::Split)]
  #[case::one(                     1.0,            true,  Expected::Target)]
  #[case::inside_child_end_band(   1.0 - 9.9e-6,   true,  Expected::Split)]
  #[case::within_rounding_of_one(  1.0 - 1e-16,    true,  Expected::Target)]
  #[case::interior_split(          0.25,           true,  Expected::Split)]
  #[case::below_half_no_split(     0.25,           false, Expected::Source)]
  #[case::above_half_no_split(     0.75,           false, Expected::Target)]
  #[trace]
  fn test_root_target_exact_placement_snaps_only_at_a_zero_distance(
    #[case] split: f64,
    #[case] split_edge: bool,
    #[case] expected: Expected,
  ) -> Result<(), Report> {
    let fixture = edge_fixture()?;
    let split = SplitFraction::new(split)?;
    let expected = expected.root_target(&fixture, split);

    let actual = fixture.target(split, Placement::Exact, split_edge)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::leaf_end(         1.0,          Expected::LeafEndSplit)]
  #[case::near_leaf_end(    1.0 - 9.9e-6, Expected::Split)]
  #[case::parent_end(       0.0,          Expected::Source)]
  #[trace]
  fn test_root_target_exact_placement_on_a_terminal_edge_is_never_the_leaf(
    #[case] split: f64,
    #[case] expected: Expected,
  ) -> Result<(), Report> {
    let fixture = terminal_edge_fixture()?;
    let split = SplitFraction::new(split)?;
    let expected = expected.root_target(&fixture, split);

    let actual = fixture.target(split, Placement::Exact, true)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_root_target_exact_placement_on_a_zero_length_edge_is_the_parent() -> Result<(), Report> {
    let fixture = zero_length_edge_fixture()?;
    let expected = RootTarget::Node(fixture.source_key);

    let actual = fixture.target(SplitFraction::MIDDLE, Placement::Exact, true)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::leaf_end_split(               1.0,            true,  Expected::LeafEndSplit)]
  #[case::inside_leaf_end_band(         1.0 - 9.9e-6,   true,  Expected::LeafEndSplit)]
  #[case::outside_leaf_end_band(        1.0 - 1.1e-5,   true,  Expected::Split)]
  #[case::interior_split(               0.75,           true,  Expected::Split)]
  #[case::parent_end(                   0.0,            true,  Expected::Source)]
  #[case::leaf_end_no_split(            1.0,            false, Expected::Source)]
  #[case::above_half_no_split(          0.75,           false, Expected::Source)]
  #[case::below_half_no_split(          0.25,           false, Expected::Source)]
  #[trace]
  fn test_root_target_on_a_terminal_edge_is_never_the_leaf(
    #[case] split: f64,
    #[case] split_edge: bool,
    #[case] expected: Expected,
  ) -> Result<(), Report> {
    let fixture = terminal_edge_fixture()?;
    let split = SplitFraction::new(split)?;
    let expected = expected.root_target(&fixture, split);

    let actual = fixture.target(split, Placement::EdgeEndBand, split_edge)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::split_at_the_root(        0.0,    true,  false)]
  #[case::inside_the_root_end_band( 9.9e-6, true,  false)]
  #[case::interior_split(           0.25,   true,  true)]
  #[case::snap_to_the_root(         0.25,   false, false)]
  #[case::snap_to_the_child(        0.75,   false, true)]
  #[trace]
  fn test_root_moves_on_a_root_edge(
    #[case] split: f64,
    #[case] split_edge: bool,
    #[case] expected: bool,
  ) -> Result<(), Report> {
    let fixture = edge_fixture()?;
    let point = EdgePoint {
      edge: fixture.edge_key,
      split: SplitFraction::new(split)?,
    };

    let actual = root_moves(
      &fixture.graph,
      &fixture.branch_lengths,
      Some(point),
      Placement::EdgeEndBand,
      split_edge,
    )?;

    assert_eq!(expected, actual);
    Ok(())
  }

  #[test]
  fn test_root_moves_without_an_edge_is_false() -> Result<(), Report> {
    let fixture = edge_fixture()?;

    let actual = root_moves(&fixture.graph, &fixture.branch_lengths, None, Placement::Exact, true)?;

    assert!(!actual);
    Ok(())
  }

  #[test]
  fn test_require_dated_new_leaves_rejects_a_new_undated_leaf() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.3,C:0.4)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let b_key = find_node_key_by_name(&graph, &names, "B").expect("B exists");
    let mut leaves_before = leaf_keys(&graph);
    leaves_before.remove(&b_key);

    let result = require_dated_new_leaves(&graph, &leaves_before, |_| false, &names);

    assert_error!(
      result,
      "Rerooting turned the internal node 'B' into a leaf, but it has no date and no observed data. Remove the node from the input tree or pass --keep-root."
    );
    Ok(())
  }

  #[test]
  fn test_require_dated_new_leaves_accepts_a_new_dated_leaf() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.3,C:0.4)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let b_key = find_node_key_by_name(&graph, &names, "B").expect("B exists");
    let mut leaves_before = leaf_keys(&graph);
    leaves_before.remove(&b_key);

    require_dated_new_leaves(&graph, &leaves_before, |key| key == b_key, &names)?;
    Ok(())
  }

  #[test]
  fn test_require_dated_new_leaves_ignores_undated_leaves_that_were_leaves_before() -> Result<(), Report> {
    let nwk_parsed = nwk_read(b"((A:0.1,B:0.2)AB:0.3,C:0.4)root;".as_slice())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;

    require_dated_new_leaves(&graph, &leaf_keys(&graph), |_| false, &names)?;
    Ok(())
  }

  mod helpers {
    use crate::reroot::find_root::EdgePoint;
    use crate::reroot::placement::{Placement, RootTarget};
    use crate::test_utils::{find_edge_key, find_node_key_by_name};
    use eyre::Report;
    use treetime_graph::key_map::EdgeMap;

    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_graph::split_fraction::SplitFraction;
    use treetime_io::nwk::nwk_read;

    #[derive(Debug, Clone, Copy)]
    pub(super) enum Expected {
      Source,
      Target,
      Split,
      LeafEndSplit,
    }

    impl Expected {
      pub(super) fn root_target(self, fixture: &EdgeFixture, split: SplitFraction) -> RootTarget {
        let edge_key = fixture.edge_key;
        match self {
          Self::Source => RootTarget::Node(fixture.source_key),
          Self::Target => RootTarget::Node(fixture.target_key),
          Self::Split => RootTarget::Split { edge_key, split },
          Self::LeafEndSplit => RootTarget::Split {
            edge_key,
            split: SplitFraction::CHILD_END,
          },
        }
      }
    }

    pub(super) struct EdgeFixture {
      pub(super) graph: Graph,
      pub(super) branch_lengths: EdgeMap<Option<f64>>,
      pub(super) edge_key: GraphEdgeKey,
      pub(super) source_key: GraphNodeKey,
      pub(super) target_key: GraphNodeKey,
    }

    impl EdgeFixture {
      pub(super) fn target(
        &self,
        split: SplitFraction,
        placement: Placement,
        split_edge: bool,
      ) -> Result<RootTarget, Report> {
        let point = EdgePoint {
          edge: self.edge_key,
          split,
        };
        RootTarget::on_edge(&self.graph, &self.branch_lengths, point, placement, split_edge)
      }
    }

    pub(super) fn edge_fixture() -> Result<EdgeFixture, Report> {
      fixture("((A:0.1,B:0.2)AB:0.3,C:0.4)root;", "root", "AB")
    }

    pub(super) fn terminal_edge_fixture() -> Result<EdgeFixture, Report> {
      fixture("((A:0.1,B:0.2)AB:0.3,C:0.4)root;", "root", "C")
    }

    pub(super) fn zero_length_edge_fixture() -> Result<EdgeFixture, Report> {
      fixture("((A:0.1,B:0.2)AB:0,C:0.4)root;", "root", "AB")
    }

    fn fixture(newick: &str, source: &str, target: &str) -> Result<EdgeFixture, Report> {
      let nwk_parsed = nwk_read(newick.as_bytes())?;
      let names = nwk_parsed.names;
      let graph = nwk_parsed.graph;
      Ok(EdgeFixture {
        edge_key: find_edge_key(&graph, &names, source, target).expect("the edge exists"),
        source_key: find_node_key_by_name(&graph, &names, source).expect("the source exists"),
        target_key: find_node_key_by_name(&graph, &names, target).expect("the target exists"),
        graph,
        branch_lengths: nwk_parsed.branch_lengths,
      })
    }
  }
}
