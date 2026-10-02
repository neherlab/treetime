#[cfg(test)]
mod tests {
  use crate::clock::clock_filter::{ClockFilterResult, clock_filter};
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_state::ClockInputs;
  use crate::clock::date_constraints::DateConstraints;
  use crate::progress::NoopProgress;
  use crate::test_utils::find_node_key_by_name;
  use crate::timetree::inference::time_inference::likely_times;
  use crate::timetree::optimization::outliers::mark_outlier_leaves;
  use eyre::Report;
  use maplit::{btreemap, btreeset};
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
    let mut date_constraints = BTreeMap::new();
    for n in graph.get_leaves() {
      let name = names.get(&n.key()).cloned().flatten();
      if let Some(name) = name {
        if let Some(&date) = dates.get(&name) {
          date_constraints.insert(n.key(), Some(Arc::new(Distribution::point(date, 1.0))));
        }
      }
    }
    DateConstraints { date_constraints }
  }

  fn clock_inputs(graph: &Graph, constraints: &DateConstraints) -> ClockInputs {
    ClockInputs::from_times(
      graph,
      &likely_times(graph, constraints, None).unwrap(),
      &BTreeMap::new(),
    )
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

    let inputs = clock_inputs(&graph, &constraints);
    let ClockFilterResult { outliers, iqd, .. } =
      clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress)?;

    assert!(outliers.is_empty(), "No outliers expected for clean data");
    assert!(iqd >= 0.0, "IQD should be non-negative");

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

    let inputs = clock_inputs(&graph, &constraints);
    let ClockFilterResult { outliers, iqd, .. } =
      clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress)?;

    assert!(iqd > 0.0, "IQD should be positive with varying dates");
    let a_key = find_node_key_by_name(&graph, &names, "A").unwrap();
    assert!(outliers.contains(&a_key), "Node A should be marked as outlier");

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

    let inputs = clock_inputs(&graph, &constraints);
    let ClockFilterResult { iqd, .. } =
      clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress)?;

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

    let inputs = clock_inputs(&graph, &constraints);
    let outliers_low_threshold = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 1.0, &NoopProgress)?
      .outliers
      .len();
    let outliers_high_threshold = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 100.0, &NoopProgress)?
      .outliers
      .len();

    assert!(
      outliers_high_threshold <= outliers_low_threshold,
      "Higher threshold should result in fewer or equal outliers"
    );

    Ok(())
  }

  #[test]
  fn test_clock_filter_mark_outlier_leaves_adds_outliers_to_leaf_flags() -> Result<(), Report> {
    let nwk_parsed = nwk_read_str("((A:0.1,B:0.1)AB:0.1,C:0.1)root;")?;
    let names = nwk_parsed.names();
    let graph = nwk_parsed.graph;
    let key = |name: &str| find_node_key_by_name(&graph, &names, name).expect("fixture node must exist");

    let outliers = btreeset! { key("B") };
    let leaf_bad_branches = btreemap! {
      key("A") => true,
      key("B") => false,
      key("C") => false,
    };

    let actual = mark_outlier_leaves(&graph, &outliers, &leaf_bad_branches);

    let expected = btreemap! {
      key("A") => true,
      key("B") => true,
      key("C") => false,
    };
    assert_eq!(expected, actual);

    Ok(())
  }
}
