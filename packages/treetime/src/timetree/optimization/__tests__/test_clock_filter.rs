#[cfg(test)]
mod tests {
  use crate::clock::clock_filter::{ClockFilterResult, clock_filter_inplace};
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_state::{ClockInputs, ClockState};
  use crate::clock::date_constraints::DateConstraints;
  use crate::timetree::optimization::clock_filter::propagate_bad_branches;
  use crate::timetree::timetree_state::TimetreeState;
  use eyre::Report;
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_distribution::Distribution;
  use treetime_graph::graph::Graph;
  use treetime_graph::node::GraphNodeKey;
  use treetime_io::nwk::nwk_read_str;

  const TREE_NEWICK: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

  fn date_constraints(
    names: &BTreeMap<GraphNodeKey, Option<String>>,
    graph: &Graph,
    dates: &BTreeMap<String, f64>,
  ) -> DateConstraints {
    let mut time_distributions = BTreeMap::new();
    for n in graph.get_leaves() {
      let name = names.get(&n.key()).cloned().flatten();
      if let Some(name) = name {
        if let Some(&date) = dates.get(&name) {
          time_distributions.insert(n.key(), Some(Arc::new(Distribution::point(date, 1.0))));
        }
      }
    }
    DateConstraints {
      time_distributions,
      ..DateConstraints::default()
    }
  }

  fn seed_clock_state(graph: &Graph, constraints: &DateConstraints) -> (ClockInputs, ClockState) {
    let date_state = TimetreeState::seed_from_values(graph, constraints);
    let inputs = ClockInputs::seed_from_times(graph, &date_state.likely_times(constraints));
    (inputs, ClockState::new(graph))
  }

  fn count_outliers(graph: &Graph, state: &ClockState) -> usize {
    graph
      .get_leaves()
      .filter(|leaf| state.node(leaf.key()).is_outlier)
      .count()
  }

  #[test]
  fn test_clock_filter_no_outliers_clean_data() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let dates = btreemap! {
      "A".to_owned() => 2010.0,
      "B".to_owned() => 2020.0,
      "C".to_owned() => 2015.0,
      "D".to_owned() => 2012.0,
    };
    let constraints = date_constraints(&names, &graph, &dates);

    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let (inputs, mut state) = seed_clock_state(&graph, &constraints);
    let ClockFilterResult { new_outliers, iqd } =
      clock_filter_inplace(&graph, &inputs, &mut state, &clock_model, &branch_lengths, 3.0)?;

    assert_eq!(count_outliers(&graph, &state), 0, "No outliers expected for clean data");
    assert!(iqd >= 0.0, "IQD should be non-negative");
    assert_eq!(new_outliers, 0, "No status changes expected");

    Ok(())
  }

  #[test]
  fn test_clock_filter_detects_outlier() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let dates = btreemap! {
      "A".to_owned() => 1900.0,
      "B".to_owned() => 2020.0,
      "C".to_owned() => 2015.0,
      "D".to_owned() => 2012.0,
    };
    let constraints = date_constraints(&names, &graph, &dates);

    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let (inputs, mut state) = seed_clock_state(&graph, &constraints);
    let ClockFilterResult { new_outliers, iqd } =
      clock_filter_inplace(&graph, &inputs, &mut state, &clock_model, &branch_lengths, 3.0)?;

    assert!(
      count_outliers(&graph, &state) >= 1,
      "At least one outlier expected for data with extreme deviation"
    );
    assert!(iqd > 0.0, "IQD should be positive with varying dates");
    assert!(new_outliers >= 1, "At least one status change expected");

    let a_is_outlier = graph.get_leaves().any(|leaf| {
      let node = leaf;
      state.node(node.key()).is_outlier
        && names
          .get(&node.key())
          .and_then(|x| x.as_deref())
          .is_some_and(|x| x == "A")
    });
    assert!(a_is_outlier, "Node A should be marked as outlier");

    Ok(())
  }

  #[test]
  fn test_clock_filter_iqd_calculation() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let dates = btreemap! {
      "A".to_owned() => 2010.0,
      "B".to_owned() => 2020.0,
      "C".to_owned() => 2015.0,
      "D".to_owned() => 2012.0,
    };
    let constraints = date_constraints(&names, &graph, &dates);

    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let (inputs, mut state) = seed_clock_state(&graph, &constraints);
    let ClockFilterResult { iqd, .. } =
      clock_filter_inplace(&graph, &inputs, &mut state, &clock_model, &branch_lengths, 3.0)?;

    assert!(iqd.is_finite(), "IQD should be a finite number");

    Ok(())
  }

  #[test]
  fn test_clock_filter_respects_threshold() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str(TREE_NEWICK)?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let branch_lengths = nwk_parsed.branch_lengths;
    let graph: Graph = graph;

    let dates = btreemap! {
      "A".to_owned() => 1980.0,
      "B".to_owned() => 2020.0,
      "C".to_owned() => 2015.0,
      "D".to_owned() => 2012.0,
    };
    let constraints = date_constraints(&names, &graph, &dates);

    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let (inputs_low, mut state_low) = seed_clock_state(&graph, &constraints);
    clock_filter_inplace(&graph, &inputs_low, &mut state_low, &clock_model, &branch_lengths, 1.0)?;
    let outliers_low_threshold = count_outliers(&graph, &state_low);

    let (inputs_high, mut state_high) = seed_clock_state(&graph, &constraints);
    clock_filter_inplace(
      &graph,
      &inputs_high,
      &mut state_high,
      &clock_model,
      &branch_lengths,
      100.0,
    )?;
    let outliers_high_threshold = count_outliers(&graph, &state_high);

    assert!(
      outliers_high_threshold <= outliers_low_threshold,
      "Higher threshold should result in fewer or equal outliers"
    );

    Ok(())
  }

  #[test]
  fn test_clock_filter_propagates_bad_branches_after_topology_change() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)AB:0.1,C:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let graph: Graph = graph;
    let mut state = TimetreeState::new(&graph);
    for node in graph.get_leaves() {
      let is_bad = names
        .get(&node.key())
        .and_then(|x| x.as_deref())
        .is_some_and(|name| name != "C");
      state.node_mut(node.key()).bad_branch = is_bad;
    }

    propagate_bad_branches(&graph, &mut state)?;

    let actual = graph
      .get_nodes()
      .map(|node| {
        let name = names[&node.key()].clone().expect("Every fixture node must be named");
        (name, state.node(node.key()).bad_branch)
      })
      .collect::<BTreeMap<_, _>>();
    let expected = btreemap! {
      "A".to_owned() => true,
      "AB".to_owned() => true,
      "B".to_owned() => true,
      "C".to_owned() => false,
      "root".to_owned() => false,
    };
    assert_eq!(expected, actual);

    Ok(())
  }
}
