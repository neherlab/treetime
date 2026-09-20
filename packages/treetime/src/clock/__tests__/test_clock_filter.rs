#![allow(
  clippy::as_conversions,
  reason = "test and benchmark code: index and expected-value casts, property-style tests over thread_rng inputs (seeding is a separate test-quality follow-up), and scratch collections"
)]

#[cfg(test)]
mod tests {
  use crate::clock::clock_filter::clock_filter_inplace;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_state::{ClockInputs, ClockState};
  use crate::o;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;
  use treetime_utils::assert_error;

  type OutlierGraphSetup = (
    Graph,
    BTreeMap<GraphNodeKey, Option<String>>,
    BTreeMap<GraphNodeKey, Option<f64>>,
    BTreeMap<GraphEdgeKey, Option<f64>>,
  );

  fn setup_outlier_graph() -> Result<OutlierGraphSetup, Report> {
    let tree = "(((A:0.1,B:0.2):0.01,(C:0.15,D:0.25):0.01):0.01,((E:0.12,F:0.18):0.01,(G:2.0,H:3.0):0.01):0.01)root;";
    let nwk_parsed = nwk_read_str(tree)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    #[rustfmt::skip]
    let dates = btreemap! {
      o!("A") => 2012.0,
      o!("B") => 2022.0,
      o!("C") => 2017.0,
      o!("D") => 2027.0,
      o!("E") => 2014.0,
      o!("F") => 2020.0,
      o!("G") => 2015.0,
      o!("H") => 2015.0,
    };

    let times = graph
      .get_leaves()
      .map(|node| {
        let name = names[&node.key()].clone().unwrap();
        (node.key(), dates.get(&name).copied())
      })
      .collect();

    Ok((graph, names, times, branch_lengths))
  }

  fn get_outlier_names(
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    graph: &Graph,
    state: &ClockState,
  ) -> Vec<String> {
    let mut result: Vec<String> = graph
      .get_leaves()
      .filter_map(|leaf| {
        let node = leaf;
        if state.node(node.key()).is_outlier {
          names.get(&node.key()).cloned().flatten()
        } else {
          None
        }
      })
      .collect();
    result.sort();
    result
  }

  #[test]
  fn test_clock_filter_positive_rate_identifies_outliers() -> Result<(), Report> {
    let (graph, names, times, branch_lengths) = setup_outlier_graph()?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let inputs = ClockInputs::seed_from_times(&graph, &times);
    let mut state = ClockState::new(&graph);
    let result = clock_filter_inplace(&graph, &inputs, &mut state, &clock_model, &branch_lengths, 3.0)?;

    assert!(result.iqd > 0.0, "IQD should be positive");
    let outliers = get_outlier_names(&names, &graph, &state);
    assert_eq!(outliers, vec![o!("G"), o!("H")]);

    Ok(())
  }

  #[test]
  fn test_clock_filter_negative_rate_identifies_same_outliers() -> Result<(), Report> {
    let (graph, names, times, branch_lengths) = setup_outlier_graph()?;
    let clock_model = ClockModel::for_testing(-0.005, 10.5);

    let inputs = ClockInputs::seed_from_times(&graph, &times);
    let mut state = ClockState::new(&graph);
    let result = clock_filter_inplace(&graph, &inputs, &mut state, &clock_model, &branch_lengths, 3.0)?;

    assert!(result.iqd > 0.0, "IQD should be positive");
    let outliers = get_outlier_names(&names, &graph, &state);
    assert_eq!(outliers, vec![o!("G"), o!("H")]);

    Ok(())
  }

  #[test]
  fn test_clock_filter_rejects_no_dated_leaves() -> Result<(), Report> {
    let (graph, times, branch_lengths) = helpers::setup_low_cardinality_graph(0)?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let inputs = ClockInputs::seed_from_times(&graph, &times);
    let mut state = ClockState::new(&graph);
    let result = clock_filter_inplace(&graph, &inputs, &mut state, &clock_model, &branch_lengths, 3.0);

    assert_error!(result, "Clock filtering requires at least one dated leaf");
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::one_dated_leaf(  1)]
  #[case::two_dated_leaves(2)]
  #[case::three_dated_leaves(3)]
  #[trace]
  fn test_clock_filter_accepts_low_cardinality_input(#[case] dated_leaf_count: usize) -> Result<(), Report> {
    let (graph, times, branch_lengths) = helpers::setup_low_cardinality_graph(dated_leaf_count)?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let inputs = ClockInputs::seed_from_times(&graph, &times);
    let mut state = ClockState::new(&graph);
    let result = clock_filter_inplace(&graph, &inputs, &mut state, &clock_model, &branch_lengths, 3.0)?;

    assert!(result.iqd.is_finite());
    Ok(())
  }

  mod helpers {
    use eyre::Report;
    use std::collections::BTreeMap;
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read_str;

    pub fn setup_low_cardinality_graph(
      dated_leaf_count: usize,
    ) -> Result<
      (
        Graph,
        BTreeMap<GraphNodeKey, Option<f64>>,
        BTreeMap<GraphEdgeKey, Option<f64>>,
      ),
      Report,
    > {
      let nwk_parsed = nwk_read_str("(A:0.1,B:0.2,C:0.3)root;")?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let branch_lengths = nwk_parsed.branch_lengths;
      let graph: Graph = graph;

      let times = graph
        .get_leaves()
        .take(dated_leaf_count)
        .enumerate()
        .map(|(index, leaf)| (leaf.key(), Some(2000.0 + index as f64)))
        .collect();

      Ok((graph, times, branch_lengths))
    }
  }
}
