#[cfg(test)]
mod tests {
  use crate::clock::clock_model::{ClockModel, ClockRegression};
  use crate::clock::clock_regression::{ClockVarianceParams, clock_regression_backward};
  use crate::clock::clock_set::ClockSet;
  use crate::clock::clock_state::{ClockInputs, ClockState};
  use crate::o;
  use crate::seq::div::{OnlyLeaves, compute_divs};
  use crate::{pretty_assert_abs_diff_eq, pretty_assert_ulps_eq};
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  fn compute_naive_rate(dates: &BTreeMap<String, f64>, div: &BTreeMap<String, f64>) -> f64 {
    let t: f64 = dates.values().sum();
    let tsq: f64 = dates.values().map(|&x| x * x).sum();
    let dt: f64 = div.iter().map(|(c, div)| div * dates[c]).sum();
    let d: f64 = div.values().sum();
    (dt * 4.0 - d * t) / (tsq * 4.0 - (t * t))
  }

  #[test]
  fn test_clock_naive_rate() -> Result<(), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    let nwk_parsed = nwk_read_str("((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;

    let graph: Graph = graph;
    let divs = compute_divs(&graph, OnlyLeaves(true), &branch_lengths, &names)?;
    let naive_rate = compute_naive_rate(&dates, &divs);

    let times = helpers::leaf_times(&names, &graph, &dates);
    let inputs = ClockInputs::seed_from_times(&graph, &times);
    let mut state = ClockState::new(&graph);
    let root_key = graph.get_exactly_one_root()?.key();

    clock_regression_backward(
      &graph,
      &inputs,
      &mut state,
      &ClockVarianceParams::default(),
      &branch_lengths,
      None,
    )?;
    let clock = ClockModel::from_regression(&ClockRegression::from_clock_set(&state.node(root_key).clock_set)?)?;
    pretty_assert_abs_diff_eq!(naive_rate, clock.clock_rate(), epsilon = 1e-10);

    let options = &ClockVarianceParams {
      variance_factor: 1.0,
      variance_offset: 0.0,
      variance_offset_leaf: 1.0,
    };

    clock_regression_backward(&graph, &inputs, &mut state, options, &branch_lengths, None)?;
    let clock = ClockModel::from_regression(&ClockRegression::from_clock_set(&state.node(root_key).clock_set)?)?;
    pretty_assert_ulps_eq!(0.007710610618916924, clock.clock_rate(), max_ulps = 4);

    Ok(())
  }

  #[test]
  fn test_clock_regression_dateless_leaf_does_not_change_root_statistics() -> Result<(), Report> {
    let dates = btreemap! {
      o!("A") => 2013.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2005.0,
    };

    let expected = helpers::root_clock_set("(A:0.1,B:0.2,C:0.2,D:0.12)root;", &dates)?;
    let actual = helpers::root_clock_set("(A:0.1,B:0.2,C:0.2,D:0.12,E:10.0)root;", &dates)?;

    assert_eq!(expected, actual);
    Ok(())
  }

  mod helpers {
    use super::*;

    pub(super) fn leaf_times(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      graph: &Graph,
      dates: &BTreeMap<String, f64>,
    ) -> BTreeMap<GraphNodeKey, Option<f64>> {
      graph
        .get_leaves()
        .map(|node| {
          let name = names[&node.key()].clone().unwrap();
          (node.key(), dates.get(&name).copied())
        })
        .collect()
    }

    pub(super) fn root_clock_set(tree: &str, dates: &BTreeMap<String, f64>) -> Result<ClockSet, Report> {
      let nwk_parsed = nwk_read_str(tree)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;
      let times = leaf_times(&names, &graph, dates);
      let inputs = ClockInputs::seed_from_times(&graph, &times);
      let mut state = ClockState::new(&graph);
      clock_regression_backward(
        &graph,
        &inputs,
        &mut state,
        &ClockVarianceParams::default(),
        &branch_lengths,
        None,
      )?;
      let root_key = graph.get_exactly_one_root()?.key();
      let clock_set = state.node(root_key).clock_set.clone();
      Ok(clock_set)
    }
  }
}
