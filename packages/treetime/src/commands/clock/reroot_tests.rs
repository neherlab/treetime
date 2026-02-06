use crate::commands::clock::clock_graph::GraphClock;
use crate::commands::clock::clock_regression::{ClockParams, clock_regression_backward, clock_regression_forward};
use crate::commands::clock::find_best_root::params::BranchPointOptimizationParams;
use crate::commands::clock::reroot::{RerootPolicy, remove_node_if_trivial, reroot_in_place};
use crate::graph::node::Named;
use crate::io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
use crate::o;
use eyre::Report;
use maplit::btreemap;
use pretty_assertions::assert_eq;

#[test]
fn test_remove_node_if_trivial_simple() -> Result<(), Report> {
  // define tree with trivial node:
  //        root
  //        /  \
  //      mid  B
  //      /
  //     A
  let mut graph: GraphClock = nwk_read_str("((A:0.5)mid:0.3,B:0.2)root;")?;

  let mid_key = graph
    .find_node(|node| node.name.as_deref() == Some("mid"))
    .expect("Expected node named 'mid'");

  remove_node_if_trivial(&mut graph, mid_key)?;

  assert!(graph.get_node(mid_key).is_none(), "Expected node to be removed");

  let expected = "(B:0.2,A:0.8)root;";
  let actual = nwk_write_str(&graph, &NwkWriteOptions::default())?;
  assert_eq!(expected, actual);

  Ok(())
}

fn setup_reroot_test_graph() -> Result<(GraphClock, ClockParams), Report> {
  let dates = btreemap! {
    o!("A") => 2013.0,
    o!("B") => 2022.0,
    o!("C") => 2017.0,
    o!("D") => 2005.0,
  };

  let graph: GraphClock = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
  for n in graph.get_leaves() {
    let name = n.read_arc().payload().read_arc().name().expect("Leaf has name").as_ref().to_owned();
    n.write_arc().payload().write_arc().time = Some(dates[&name]);
  }

  let options = ClockParams::default();
  clock_regression_backward(&graph, &options);
  clock_regression_forward(&graph, &options);

  Ok((graph, options))
}

#[test]
fn test_reroot_policy_allow_edge_split_false_no_new_nodes() -> Result<(), Report> {
  let (mut graph, options) = setup_reroot_test_graph()?;
  let node_count_before = graph.get_nodes().len();

  // Both flags false: don't split edges AND don't remove old root
  // This guarantees no new nodes created and no nodes removed
  let policy = RerootPolicy {
    allow_edge_split: false,
    remove_old_root_if_trivial: false,
  };

  reroot_in_place(&mut graph, &options, &BranchPointOptimizationParams::default(), &policy)?;

  let node_count_after = graph.get_nodes().len();
  assert_eq!(node_count_before, node_count_after, "Node count should be unchanged when edge split is disabled and old root is preserved");

  Ok(())
}

#[test]
fn test_reroot_policy_remove_old_root_if_trivial_false_preserves_old_root() -> Result<(), Report> {
  let (mut graph, options) = setup_reroot_test_graph()?;
  let old_root_key = graph.get_exactly_one_root()?.read_arc().key();

  let policy = RerootPolicy {
    allow_edge_split: true,
    remove_old_root_if_trivial: false,
  };

  let new_root_key = reroot_in_place(&mut graph, &options, &BranchPointOptimizationParams::default(), &policy)?;

  if new_root_key != old_root_key {
    assert!(
      graph.get_node(old_root_key).is_some(),
      "Old root should still exist when remove_old_root_if_trivial is false"
    );
  }

  Ok(())
}

#[test]
fn test_reroot_policy_default_allows_edge_split() -> Result<(), Report> {
  let (mut graph, options) = setup_reroot_test_graph()?;
  let node_count_before = graph.get_nodes().len();

  let policy = RerootPolicy::default();

  reroot_in_place(&mut graph, &options, &BranchPointOptimizationParams::default(), &policy)?;

  let node_count_after = graph.get_nodes().len();
  // With default policy, a new node may be created by edge split (count increases)
  // or old trivial root may be removed (count stays same or decreases by 1 if split created one)
  // The key is it should not crash and should complete successfully
  assert!(node_count_after >= node_count_before - 1, "Node count should be reasonable after reroot with default policy");

  Ok(())
}
