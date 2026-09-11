#[cfg(test)]
mod tests {
  use crate::clock::clock_graph::GraphClock;
  use crate::clock::clock_state::{ClockEdgeState, ClockNodeState, ClockState};
  use crate::payload::clock_set::ClockSet;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  // Oracle: the value inputs passed to `seed_from_values` plus its documented contract -- the seeded
  // date comes from `times` (`None` for a missing key), and `div`/`is_outlier`/`bad_branch`/`clock_set`
  // and every edge start at their defaults.
  #[test]
  fn test_clock_state_seed_from_values_sources_times_and_defaults_the_rest() -> Result<(), Report> {
    let graph: GraphClock = nwk_read_str("(A:0.1,B:0.2)root;")?;
    let key_of = helpers::key_by_name(&graph);
    let (a, b, root) = (key_of["A"], key_of["B"], key_of["root"]);

    // Seed A with a date, leave B explicitly dateless, and omit root entirely (a missing key).
    let times = btreemap! {
      a => Some(2000.0),
      b => None,
    };

    let state = ClockState::seed_from_values(&graph, &times);

    let expected_nodes = btreemap! {
      a    => helpers::node_with_time(Some(2000.0)),
      b    => helpers::node_with_time(None),
      root => helpers::node_with_time(None),
    };
    assert_eq!(expected_nodes, state.nodes);

    let expected_edges: BTreeMap<_, _> = graph
      .get_edges()
      .iter()
      .map(|edge| (edge.read_arc().key(), ClockEdgeState::default()))
      .collect();
    assert_eq!(expected_edges, state.edges);

    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn key_by_name(graph: &GraphClock) -> BTreeMap<String, GraphNodeKey> {
      graph
        .get_nodes()
        .iter()
        .map(|node| {
          let node = node.read_arc();
          let name = node.payload().read_arc().name.clone().expect("node has a name");
          (name, node.key())
        })
        .collect()
    }

    pub(super) fn node_with_time(time: Option<f64>) -> ClockNodeState {
      ClockNodeState {
        clock_set: ClockSet::default(),
        div: 0.0,
        time,
        bad_branch: false,
        is_outlier: false,
      }
    }
  }
}
