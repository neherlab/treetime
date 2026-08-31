#[cfg(test)]
mod tests {
  use crate::partition::timetree::partition::GraphTimetree;
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use crate::timetree::optimization::polytomy::apply::{ChildRef, apply_plan};
  use crate::timetree::optimization::polytomy::sweep::{Merger, SubtreePlan};
  use eyre::Report;
  use pretty_assertions::assert_eq;
  use treetime_graph::edge::{GraphEdgeKey, HasBranchLength};
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::io::json::{JsonPretty, json_write_str};
  use treetime_utils::{assert_error, make_report, pretty_assert_abs_diff_eq};

  /// `((A,B,C)P)root` with times set and each child edge carrying a distinct mutation length.
  fn polytomy_graph() -> Result<(GraphTimetree, GraphNodeKey, Vec<ChildRef>), Report> {
    let graph: GraphTimetree = nwk_read_str("((A:0.1,B:0.2,C:0.15)P:0.05)root;")?;

    let times = [
      ("A", 2020.0),
      ("B", 2018.0),
      ("C", 2016.0),
      ("P", 2000.0),
      ("root", 1990.0),
    ];
    for (name, time) in times {
      let key = find_node_key_by_name(&graph, name).ok_or_else(|| make_report!("{name} not found"))?;
      graph
        .get_node(key)
        .expect("Node must exist")
        .write_arc()
        .payload()
        .write_arc()
        .time = Some(time);
    }

    let parent_key = find_node_key_by_name(&graph, "P").ok_or_else(|| make_report!("P not found"))?;

    let mut children = Vec::new();
    for (name, time, mutation_length) in [("A", 2020.0, 0.3), ("B", 2018.0, 0.4), ("C", 2016.0, 0.5)] {
      let node_key = find_node_key_by_name(&graph, name).ok_or_else(|| make_report!("{name} not found"))?;
      let edge_key = find_edge_key(&graph, "P", name).ok_or_else(|| make_report!("P->{name} not found"))?;
      let edge = graph.get_edge(edge_key).expect("Edge must exist");
      edge
        .write_arc()
        .payload()
        .write_arc()
        .set_branch_length(Some(mutation_length));
      children.push(ChildRef {
        node_key,
        edge_key,
        time,
      });
    }

    Ok((graph, parent_key, children))
  }

  fn child_edge_of(graph: &GraphTimetree, node_key: GraphNodeKey) -> GraphEdgeKey {
    let node = graph.get_node(node_key).expect("Node must exist");
    let node = node.read_arc();
    assert_eq!(node.inbound().len(), 1, "a tree node has exactly one parent edge");
    node.inbound()[0]
  }

  fn parent_of(graph: &GraphTimetree, node_key: GraphNodeKey) -> GraphNodeKey {
    let edge_key = child_edge_of(graph, node_key);
    graph.get_edge(edge_key).expect("Edge must exist").read_arc().source()
  }

  #[test]
  fn test_apply_plan_reparents_children_keeping_edge_key_and_mutation_length() -> Result<(), Report> {
    let (mut graph, parent_key, children) = polytomy_graph()?;
    let original_edges: Vec<GraphEdgeKey> = children.iter().map(|child| child.edge_key).collect();

    // Merge A and B (lineages 0 and 1) at 2010; C stays a direct child of P.
    let plan = SubtreePlan {
      mergers: vec![Merger {
        time: 2010.0,
        left: 0,
        right: 1,
      }],
      roots: vec![2, 3],
    };

    let created = apply_plan(&mut graph, parent_key, 2000.0, &children, &plan)?;
    assert_eq!(created, 1);

    for (index, child) in children.iter().enumerate() {
      assert_eq!(
        child_edge_of(&graph, child.node_key),
        original_edges[index],
        "reparenting must preserve the edge key so partition state keyed by it survives"
      );
    }

    for (child, expected) in children.iter().zip([0.3, 0.4, 0.5]) {
      let edge = graph.get_edge(child.edge_key).expect("Edge must exist");
      let mutation_length = edge.read_arc().payload().read_arc().branch_length();
      pretty_assert_abs_diff_eq!(
        mutation_length.expect("observed mutation length must survive reparenting"),
        expected,
        epsilon = 1e-12
      );
    }

    Ok(())
  }

  #[test]
  fn test_apply_plan_builds_the_planned_topology() -> Result<(), Report> {
    let (mut graph, parent_key, children) = polytomy_graph()?;

    let plan = SubtreePlan {
      mergers: vec![Merger {
        time: 2010.0,
        left: 0,
        right: 1,
      }],
      roots: vec![2, 3],
    };
    apply_plan(&mut graph, parent_key, 2000.0, &children, &plan)?;

    let merger_key = parent_of(&graph, children[0].node_key);
    assert_eq!(
      parent_of(&graph, children[1].node_key),
      merger_key,
      "both merged children must hang off the same new node"
    );
    assert_ne!(
      merger_key, parent_key,
      "the merged pair must not stay under the polytomy"
    );
    assert_eq!(
      parent_of(&graph, children[2].node_key),
      parent_key,
      "an unmerged child stays a direct child of the polytomy"
    );
    assert_eq!(parent_of(&graph, merger_key), parent_key);

    let parent_degree = graph
      .get_node(parent_key)
      .expect("Node must exist")
      .read_arc()
      .degree_out();
    assert_eq!(
      parent_degree, 2,
      "a 3-way polytomy with one merger becomes a bifurcation"
    );

    pretty_assert_abs_diff_eq!(
      graph
        .get_node(merger_key)
        .expect("Node must exist")
        .read_arc()
        .payload()
        .read_arc()
        .time
        .expect("new node must be dated"),
      2010.0,
      epsilon = 1e-12
    );

    Ok(())
  }

  #[test]
  fn test_apply_plan_sets_time_lengths_from_the_new_parent() -> Result<(), Report> {
    let (mut graph, parent_key, children) = polytomy_graph()?;

    let plan = SubtreePlan {
      mergers: vec![Merger {
        time: 2010.0,
        left: 0,
        right: 1,
      }],
      roots: vec![2, 3],
    };
    apply_plan(&mut graph, parent_key, 2000.0, &children, &plan)?;

    let time_length_of = |edge_key: GraphEdgeKey| -> f64 {
      graph
        .get_edge(edge_key)
        .expect("Edge must exist")
        .read_arc()
        .payload()
        .read_arc()
        .time_length
        .expect("time_length must be set")
    };

    // A and B now hang off the merger at 2010.
    pretty_assert_abs_diff_eq!(time_length_of(children[0].edge_key), 10.0, epsilon = 1e-12);
    pretty_assert_abs_diff_eq!(time_length_of(children[1].edge_key), 8.0, epsilon = 1e-12);
    // C still hangs off P at 2000.
    pretty_assert_abs_diff_eq!(time_length_of(children[2].edge_key), 16.0, epsilon = 1e-12);

    let merger_key = parent_of(&graph, children[0].node_key);
    let merger_edge = child_edge_of(&graph, merger_key);
    pretty_assert_abs_diff_eq!(time_length_of(merger_edge), 10.0, epsilon = 1e-12);

    Ok(())
  }

  #[test]
  fn test_apply_plan_gives_merger_edges_zero_mutation_length() -> Result<(), Report> {
    let (mut graph, parent_key, children) = polytomy_graph()?;

    let plan = SubtreePlan {
      mergers: vec![Merger {
        time: 2010.0,
        left: 0,
        right: 1,
      }],
      roots: vec![2, 3],
    };
    apply_plan(&mut graph, parent_key, 2000.0, &children, &plan)?;

    let merger_key = parent_of(&graph, children[0].node_key);
    let merger_edge = child_edge_of(&graph, merger_key);
    let mutation_length = graph
      .get_edge(merger_edge)
      .expect("Edge must exist")
      .read_arc()
      .payload()
      .read_arc()
      .branch_length();

    pretty_assert_abs_diff_eq!(
      mutation_length.expect("a new branch must have an explicit mutation length"),
      0.0,
      epsilon = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_apply_plan_nests_mergers_that_reference_earlier_mergers() -> Result<(), Report> {
    let (mut graph, parent_key, children) = polytomy_graph()?;

    // ((A,B) at 2012, C) at 2005 -- the second merger consumes the first.
    let plan = SubtreePlan {
      mergers: vec![
        Merger {
          time: 2012.0,
          left: 0,
          right: 1,
        },
        Merger {
          time: 2005.0,
          left: 2,
          right: 3,
        },
      ],
      roots: vec![4],
    };

    let created = apply_plan(&mut graph, parent_key, 2000.0, &children, &plan)?;
    assert_eq!(created, 2);

    let inner = parent_of(&graph, children[0].node_key);
    let outer = parent_of(&graph, inner);
    assert_eq!(parent_of(&graph, children[1].node_key), inner);
    assert_eq!(parent_of(&graph, children[2].node_key), outer);
    assert_eq!(parent_of(&graph, outer), parent_key);

    let parent_degree = graph
      .get_node(parent_key)
      .expect("Node must exist")
      .read_arc()
      .degree_out();
    assert_eq!(
      parent_degree, 1,
      "a fully consumed polytomy leaves one child for cleanup"
    );

    Ok(())
  }

  #[test]
  fn test_apply_plan_rejects_a_merger_referencing_a_later_merger() -> Result<(), Report> {
    let (mut graph, parent_key, children) = polytomy_graph()?;

    // Lineage 4 is the second merger's own node, which does not exist yet.
    let plan = SubtreePlan {
      mergers: vec![
        Merger {
          time: 2012.0,
          left: 0,
          right: 4,
        },
        Merger {
          time: 2005.0,
          left: 1,
          right: 2,
        },
      ],
      roots: vec![3],
    };
    let before = json_write_str(&graph, JsonPretty(false))?;

    assert_error!(
      apply_plan(&mut graph, parent_key, 2000.0, &children, &plan),
      "Polytomy plan merger 0 referenced lineage 4 before it was created. This is an internal error. Please report it to developers."
    );
    let after = json_write_str(&graph, JsonPretty(false))?;
    assert_eq!(before, after, "plan validation must finish before graph mutation");
    Ok(())
  }
}
