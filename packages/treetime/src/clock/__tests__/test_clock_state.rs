#[cfg(test)]
mod tests {
  use crate::clock::clock_state::{ClockEdgeInput, ClockInputs, ClockNodeInput};
  use crate::test_utils::{find_edge_key, find_node_key_by_name};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_io::nwk::nwk_read_str;

  #[test]
  fn test_clock_inputs_from_times_sources_times_and_edge_inputs_and_defaults_the_rest() -> Result<(), Report> {
    let parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let names = parsed.names();
    let graph = parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).unwrap();
    let (a, b, root) = (key("A"), key("B"), key("root"));
    let root_a = find_edge_key(&graph, &names, "root", "A").unwrap();
    let root_b = find_edge_key(&graph, &names, "root", "B").unwrap();

    let times = btreemap! {
      a => Some(2000.0),
      b => None,
    };
    let edge_inputs = btreemap! {
      root_a => (Some(1.5), 2.0),
    };

    let inputs = ClockInputs::from_times(&graph, &times, &edge_inputs);

    let node = |time| ClockNodeInput {
      time,
      bad_branch: false,
    };
    let expected_nodes = btreemap! {
      a    => node(Some(2000.0)),
      b    => node(None),
      root => node(None),
    };
    assert_eq!(expected_nodes, inputs.nodes);

    let expected_edges: BTreeMap<_, _> = btreemap! {
      root_a => ClockEdgeInput { time_length: Some(1.5), gamma: 2.0 },
      root_b => ClockEdgeInput::default(),
    };
    assert_eq!(expected_edges, inputs.edges);

    Ok(())
  }
}
