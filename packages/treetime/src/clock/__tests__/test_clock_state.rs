#[cfg(test)]
mod tests {
  use crate::clock::clock_state::{ClockEdgeInput, ClockInputs, ClockNodeInput};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  // Oracle: the value inputs passed to `seed_from_times` plus its documented contract -- the seeded
  // date comes from `times` (`None` for a missing key), and `bad_branch` and every edge start at
  // their defaults.
  #[test]
  fn test_clock_state_seed_from_times_sources_times_and_defaults_the_rest() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let key_of = helpers::key_by_name(&names, &graph);
    let (a, b, root) = (key_of["A"], key_of["B"], key_of["root"]);

    // Seed A with a date, leave B explicitly dateless, and omit root entirely (a missing key).
    let times = btreemap! {
      a => Some(2000.0),
      b => None,
    };

    let inputs = ClockInputs::seed_from_times(&graph, &times);

    let expected_nodes = btreemap! {
      a    => helpers::node_with_time(Some(2000.0)),
      b    => helpers::node_with_time(None),
      root => helpers::node_with_time(None),
    };
    assert_eq!(expected_nodes, inputs.nodes);

    let expected_edges: BTreeMap<_, _> = graph
      .get_edges()
      .map(|edge| (edge.key(), ClockEdgeInput::default()))
      .collect();
    assert_eq!(expected_edges, inputs.edges);

    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn key_by_name(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
    ) -> BTreeMap<String, GraphNodeKey> {
      graph
        .get_nodes()
        .map(|node| {
          let name = names.get(&node.key()).cloned().flatten().expect("node has a name");
          (name, node.key())
        })
        .collect()
    }

    pub(super) fn node_with_time(time: Option<f64>) -> ClockNodeInput {
      ClockNodeInput {
        time,
        bad_branch: false,
      }
    }
  }
}
