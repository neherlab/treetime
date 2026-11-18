#[cfg(test)]
mod tests {
  use crate::graph::Graph;
  use crate::key_map::{EdgeMap, NodeMap};
  use crate::merge_by_name::{MergedGraph, NamedTree, merge_trees_by_name};
  use crate::node::GraphNodeKey;
  use crate::node_names::NodeNames;
  use crate::pass_through::PassThroughReport;
  use crate::schedule::GraphSchedule;
  use crate::scope::SiteScope;
  use eyre::Report;
  use indexmap::IndexMap;
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeSet;
  use treetime_utils::assert_error;

  struct TestTree {
    label: String,
    graph: Graph,
    names: NodeNames,
    unnamed: BTreeSet<GraphNodeKey>,
    lengths: EdgeMap<Option<f64>>,
  }

  impl TestTree {
    fn new(label: &str, edges: &[(&str, &str, f64)]) -> Result<Self, Report> {
      let mut graph = Graph::new();
      let mut keys: IndexMap<&str, GraphNodeKey> = IndexMap::new();
      let mut input: NodeMap<Option<String>> = NodeMap::new();
      let mut lengths = EdgeMap::new();
      let mut unnamed_count = 0;
      for &(parent, child, length) in edges {
        let [parent_key, child_key] = [parent, child].map(|name| {
          if name.is_empty() {
            unnamed_count += 1;
            let key = graph.add_node();
            input.insert(key, None);
            return key;
          }
          *keys.entry(name).or_insert_with(|| {
            let key = graph.add_node();
            input.insert(key, Some(name.to_owned()));
            key
          })
        });
        let edge = graph.add_edge(parent_key, child_key)?;
        lengths.insert(edge, Some(length));
      }
      graph.build()?;
      let parsed = NodeNames::from_input(input, &graph)?;
      assert_eq!(unnamed_count, parsed.unnamed_nodes.len());
      Ok(Self {
        label: label.to_owned(),
        graph,
        names: parsed.names,
        unnamed: parsed.unnamed_nodes,
        lengths,
      })
    }

    fn named(&self) -> NamedTree<'_> {
      NamedTree {
        label: &self.label,
        graph: &self.graph,
        names: &self.names,
        unnamed_nodes: &self.unnamed,
        branch_lengths: &self.lengths,
      }
    }
  }

  fn merge(trees: &[TestTree]) -> Result<MergedGraph, Report> {
    merge_trees_by_name(&trees.iter().map(TestTree::named).collect_vec())
  }

  fn edges(merged: &MergedGraph) -> Vec<(String, String, Option<f64>, Vec<usize>)> {
    merged
      .graph
      .get_edges()
      .map(|edge| {
        let segments = (0..merged.scopes.segment_count())
          .filter(|&segment| merged.scopes.carries(edge.key(), segment))
          .collect_vec();
        (
          merged.names[edge.source()].to_owned(),
          merged.names[edge.target()].to_owned(),
          merged.branch_lengths[&edge.key()],
          segments,
        )
      })
      .collect()
  }

  fn names(merged: &MergedGraph, keys: &[GraphNodeKey]) -> Vec<String> {
    keys.iter().map(|&key| merged.names[key].to_owned()).collect()
  }

  fn reassortment() -> Result<Vec<TestTree>, Report> {
    Ok(vec![
      TestTree::new(
        "seg1",
        &[
          ("R", "X", 0.1),
          ("R", "Y", 0.2),
          ("X", "A", 0.3),
          ("X", "C", 0.1),
          ("C", "c1", 0.1),
          ("C", "c2", 0.1),
          ("Y", "B1", 0.2),
          ("Y", "B2", 0.2),
        ],
      )?,
      TestTree::new(
        "seg2",
        &[
          ("R", "X", 0.3),
          ("R", "Y", 0.2),
          ("X", "A", 0.3),
          ("Y", "C", 0.0),
          ("C", "c1", 0.1),
          ("C", "c2", 0.1),
          ("Y", "B1", 0.2),
          ("Y", "B2", 0.2),
        ],
      )?,
    ])
  }

  #[test]
  fn test_merge_trees_by_name_shares_nodes_and_edges_and_records_segments() -> Result<(), Report> {
    let merged = merge(&reassortment()?)?;
    let s = |name: &str| name.to_owned();
    assert_eq!(
      vec![
        (s("R"), s("X"), Some(0.2), vec![0, 1]),
        (s("X"), s("A"), Some(0.3), vec![0, 1]),
        (s("X"), s("C"), Some(0.1), vec![0]),
        (s("C"), s("c1"), Some(0.1), vec![0, 1]),
        (s("C"), s("c2"), Some(0.1), vec![0, 1]),
        (s("R"), s("Y"), Some(0.2), vec![0, 1]),
        (s("Y"), s("B1"), Some(0.2), vec![0, 1]),
        (s("Y"), s("B2"), Some(0.2), vec![0, 1]),
        (s("Y"), s("C"), Some(0.0), vec![1]),
      ],
      edges(&merged)
    );
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_gives_a_reassortant_node_one_parent_per_segment() -> Result<(), Report> {
    let merged = merge(&reassortment()?)?;
    let schedule = GraphSchedule::new(&merged.graph)?;
    let c = merged
      .graph
      .get_nodes()
      .find(|node| &merged.names[node.key()] == "C")
      .map(|node| node.key());
    let parents = schedule.parents(c.expect("C exists"));
    let by_segment = (0..2)
      .map(|segment| {
        parents
          .iter()
          .filter(|(_, edge)| merged.scopes.get(*edge, segment) == &SiteScope::Full)
          .map(|(parent, _)| merged.names[*parent].to_owned())
          .collect_vec()
      })
      .collect_vec();
    assert_eq!(vec![vec!["X".to_owned()], vec!["Y".to_owned()]], by_segment);
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_keeps_several_roots() -> Result<(), Report> {
    let trees = vec![
      TestTree::new("seg1", &[("R", "X", 0.1), ("X", "A", 0.1), ("X", "B", 0.1)])?,
      TestTree::new("seg2", &[("R2", "X", 0.1), ("X", "A", 0.1), ("X", "B", 0.1)])?,
    ];
    let merged = merge(&trees)?;
    let schedule = GraphSchedule::new(&merged.graph)?;
    assert_eq!(vec!["R", "R2"], names(&merged, schedule.roots()));
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_averages_the_lengths_of_a_shared_edge_over_the_trees_that_give_one() -> Result<(), Report>
  {
    let mut without_lengths = TestTree::new("seg3", &[("R", "A", 0.0), ("R", "B", 0.0)])?;
    for length in without_lengths.lengths.values_mut() {
      *length = None;
    }
    let trees = vec![
      TestTree::new("seg1", &[("R", "A", 0.1), ("R", "B", 0.3)])?,
      TestTree::new("seg2", &[("R", "A", 0.2), ("R", "B", 0.1)])?,
      without_lengths,
    ];
    let lengths = edges(&merge(&trees)?)
      .into_iter()
      .map(|(_, child, length, _)| (child, length.map(|length| (length * 1e12).round() / 1e12)))
      .collect_vec();
    assert_eq!(vec![("A".to_owned(), Some(0.15)), ("B".to_owned(), Some(0.2))], lengths);
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_keeps_no_length_when_no_tree_gives_one() -> Result<(), Report> {
    let mut tree = TestTree::new("seg1", &[("R", "A", 0.0), ("R", "B", 0.0)])?;
    for length in tree.lengths.values_mut() {
      *length = None;
    }
    let lengths = edges(&merge(&[tree])?).into_iter().map(|edge| edge.2).collect_vec();
    assert_eq!(vec![None, None], lengths);
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_rejects_conflicting_ancestry() -> Result<(), Report> {
    let trees = vec![
      TestTree::new(
        "seg1",
        &[("R", "X", 0.1), ("X", "Y", 0.1), ("Y", "A", 0.1), ("R", "B", 0.1)],
      )?,
      TestTree::new(
        "seg2",
        &[("R", "Y", 0.1), ("Y", "X", 0.1), ("X", "A", 0.1), ("R", "B", 0.1)],
      )?,
    ];
    assert_error!(
      merge(&trees),
      "The trees order their nodes in conflicting ways: merged by name they form the cycle 'Y' -> 'X' -> 'Y', where each node is a parent of the next"
    );
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_rejects_unnamed_nodes() -> Result<(), Report> {
    let trees = vec![TestTree::new("seg1", &[("R", "", 0.1), ("R", "B", 0.1)])?];
    assert_error!(
      merge(&trees),
      "The tree 'seg1' has 1 nodes without a name. Trees merge by node name, so every node needs one"
    );
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_rejects_a_name_repeated_within_a_tree() -> Result<(), Report> {
    let mut graph = Graph::new();
    let [r, a, b] = [(); 3].map(|()| graph.add_node());
    let ra = graph.add_edge(r, a)?;
    let rb = graph.add_edge(r, b)?;
    graph.build()?;
    let input: NodeMap<Option<String>> = [(r, "R"), (a, "A"), (b, "A")]
      .into_iter()
      .map(|(key, name)| (key, Some(name.to_owned())))
      .collect();
    let names = NodeNames::from_input(input, &graph)?.names;
    let lengths: EdgeMap<Option<f64>> = [(ra, Some(0.1)), (rb, Some(0.1))].into_iter().collect();
    let tree = NamedTree {
      label: "seg1",
      graph: &graph,
      names: &names,
      unnamed_nodes: &BTreeSet::new(),
      branch_lengths: &lengths,
    };
    assert_error!(
      merge_trees_by_name(&[tree]),
      "The tree 'seg1' gives the same name to more than one node: 'A'. Trees merge by node name, so names must be unique within a tree"
    );
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_routes_each_segment_through_the_nodes_of_the_other_and_splits_its_length()
  -> Result<(), Report> {
    let trees = vec![
      TestTree::new(
        "seg1",
        &[("R", "X", 0.1), ("X", "A", 0.3), ("X", "B", 0.3), ("R", "C", 0.4)],
      )?,
      TestTree::new(
        "seg2",
        &[("R", "A", 0.8), ("R", "BC", 0.1), ("BC", "B", 0.3), ("BC", "C", 0.3)],
      )?,
    ];
    let merged = merge(&trees)?;
    let s = |name: &str| name.to_owned();
    let rounded = edges(&merged)
      .into_iter()
      .map(|(parent, child, length, segments)| {
        (
          parent,
          child,
          length.map(|length| (length * 1e9).round() / 1e9),
          segments,
        )
      })
      .collect_vec();
    assert_eq!(
      (
        PassThroughReport {
          inserted_nodes: 2,
          ambiguous_edges: 0,
          conflicting_edges: 0,
        },
        vec![
          (s("R"), s("X"), Some(0.15), vec![0, 1]),
          (s("X"), s("A"), Some(0.45), vec![0, 1]),
          (s("X"), s("B"), Some(0.3), vec![0]),
          (s("R"), s("BC"), Some(0.1), vec![0, 1]),
          (s("BC"), s("C"), Some(0.3), vec![0, 1]),
          (s("BC"), s("B"), Some(0.3), vec![1]),
        ],
      ),
      (merged.pass_through, rounded)
    );
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_leaves_an_edge_whose_other_segments_disagree_on_the_path() -> Result<(), Report> {
    let trees = vec![
      TestTree::new("seg1", &[("R", "W", 0.1), ("W", "A", 0.1), ("W", "B", 0.1)])?,
      TestTree::new("seg2", &[("R", "V", 0.1), ("V", "A", 0.1), ("V", "B", 0.1)])?,
      TestTree::new("seg3", &[("R", "A", 0.2), ("R", "B", 0.2)])?,
    ];
    let merged = merge(&trees)?;
    assert_eq!(
      (
        PassThroughReport {
          inserted_nodes: 0,
          ambiguous_edges: 2,
          conflicting_edges: 0,
        },
        vec![2],
      ),
      (merged.pass_through, segments_of(&merged, "R", "A"))
    );
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_leaves_an_edge_whose_path_has_a_node_the_segment_holds_elsewhere() -> Result<(), Report> {
    let trees = vec![
      TestTree::new(
        "seg1",
        &[("R", "W", 0.1), ("W", "A", 0.1), ("W", "C", 0.1), ("R", "B", 0.2)],
      )?,
      TestTree::new(
        "seg2",
        &[("R", "A", 0.2), ("R", "W", 0.1), ("W", "B", 0.1), ("W", "C", 0.1)],
      )?,
    ];
    let merged = merge(&trees)?;
    assert_eq!(
      (
        PassThroughReport {
          inserted_nodes: 0,
          ambiguous_edges: 0,
          conflicting_edges: 2,
        },
        vec![1],
        vec![0],
      ),
      (
        merged.pass_through,
        segments_of(&merged, "R", "A"),
        segments_of(&merged, "R", "B")
      )
    );
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_routes_every_sibling_of_a_polytomy_through_the_node_another_segment_puts_above_them()
  -> Result<(), Report> {
    let trees = vec![
      TestTree::new(
        "seg1",
        &[
          ("X", "W", 0.1),
          ("W", "Y1", 0.1),
          ("W", "Y2", 0.1),
          ("W", "Y3", 0.1),
          ("X", "Z", 0.2),
        ],
      )?,
      TestTree::new(
        "seg2",
        &[("X", "Y1", 0.2), ("X", "Y2", 0.2), ("X", "Y3", 0.2), ("X", "Z", 0.2)],
      )?,
    ];
    let merged = merge(&trees)?;
    let s = |name: &str| name.to_owned();
    assert_eq!(
      (
        PassThroughReport {
          inserted_nodes: 1,
          ambiguous_edges: 0,
          conflicting_edges: 0,
        },
        vec![
          (s("X"), s("W"), Some(0.1), vec![0, 1]),
          (s("W"), s("Y1"), Some(0.1), vec![0, 1]),
          (s("W"), s("Y2"), Some(0.1), vec![0, 1]),
          (s("W"), s("Y3"), Some(0.1), vec![0, 1]),
          (s("X"), s("Z"), Some(0.2), vec![0, 1]),
        ],
      ),
      (merged.pass_through, edges(&merged))
    );
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_routes_an_edge_through_a_node_the_segment_passes_through_with_one_child()
  -> Result<(), Report> {
    let trees = vec![
      TestTree::new("seg1", &[("X", "W", 0.1), ("W", "c", 0.1), ("W", "Y", 0.1)])?,
      TestTree::new("seg2", &[("X", "W", 0.1), ("W", "c", 0.1), ("X", "Y", 0.2)])?,
    ];
    let merged = merge(&trees)?;
    assert_eq!(
      (0, vec![0, 1], Vec::<usize>::new()),
      (
        merged.pass_through.conflicting_edges,
        segments_of(&merged, "W", "Y"),
        segments_of(&merged, "X", "Y"),
      )
    );
    Ok(())
  }

  #[test]
  fn test_merge_trees_by_name_refines_paths_until_every_segment_shares_the_longest() -> Result<(), Report> {
    let trees = vec![
      TestTree::new(
        "seg1",
        &[
          ("X", "W", 0.1),
          ("W", "U", 0.1),
          ("U", "Y", 0.1),
          ("W", "a", 0.1),
          ("U", "b", 0.1),
        ],
      )?,
      TestTree::new(
        "seg2",
        &[("X", "W", 0.1), ("W", "Y", 0.2), ("W", "a", 0.1), ("X", "b", 0.3)],
      )?,
      TestTree::new("seg3", &[("X", "Y", 0.3), ("X", "a", 0.2), ("X", "b", 0.3)])?,
    ];
    let merged = merge(&trees)?;
    assert_eq!(
      (vec![0, 1, 2], vec![0, 1, 2], vec![0, 1, 2], 3),
      (
        segments_of(&merged, "X", "W"),
        segments_of(&merged, "W", "U"),
        segments_of(&merged, "U", "Y"),
        merged.pass_through.inserted_nodes,
      )
    );
    Ok(())
  }

  fn segments_of(merged: &MergedGraph, parent: &str, child: &str) -> Vec<usize> {
    edges(merged)
      .into_iter()
      .find(|(p, c, _, _)| p == parent && c == child)
      .map(|(_, _, _, segments)| segments)
      .unwrap_or_default()
  }
}
