use crate::commands::clock::clock_graph::GraphClock;
use crate::commands::clock::reroot::remove_node_if_trivial;
use crate::io::nwk::{NwkWriteOptions, nwk_read_str, nwk_write_str};
use eyre::Report;
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
