#[cfg(test)]
mod tests {
  use crate::clock::clock_filter::clock_filter;
  use crate::clock::clock_model::ClockModel;
  use crate::clock::clock_state::ClockInputs;
  use crate::o;
  use crate::progress::NoopProgress;
  use eyre::Report;
  use helpers::{
    IQD_TREE, OUTLIER_TREE, divergences_by_name, iqd_dates, outlier_dates, outlier_names, setup_graph,
    setup_low_cardinality_graph,
  };
  use maplit::btreemap;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use std::collections::BTreeMap;
  use treetime_utils::assert_error;
  use treetime_utils::{pretty_assert_abs_diff_eq, pretty_assert_map_abs_diff_eq};

  #[test]
  fn test_clock_filter_positive_rate_identifies_outliers() -> Result<(), Report> {
    let (graph, names, times, branch_lengths) = setup_graph(OUTLIER_TREE, &outlier_dates())?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let result = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress)?;

    assert_eq!(vec![o!("G"), o!("H")], outlier_names(&names, &result.outliers));
    pretty_assert_map_abs_diff_eq!(
      helpers::outlier_tree_divergences(),
      divergences_by_name(&names, &result.divergences),
      epsilon = 1e-12
    );
    Ok(())
  }

  #[test]
  fn test_clock_filter_negative_rate_identifies_same_outliers() -> Result<(), Report> {
    let (graph, names, times, branch_lengths) = setup_graph(OUTLIER_TREE, &outlier_dates())?;
    let clock_model = ClockModel::for_testing(-0.005, 10.5);

    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let result = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress)?;

    assert_eq!(vec![o!("G"), o!("H")], outlier_names(&names, &result.outliers));
    pretty_assert_map_abs_diff_eq!(
      helpers::outlier_tree_divergences(),
      divergences_by_name(&names, &result.divergences),
      epsilon = 1e-12
    );
    Ok(())
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::boundary_is_not_an_outlier( 3.0, vec![o!("I")])]
  #[case::lower_threshold(            0.5, vec![o!("A"), o!("F"), o!("G"), o!("H"), o!("I")])]
  #[trace]
  fn test_clock_filter_flags_deviations_beyond_the_threshold_times_the_interquartile_distance(
    #[case] threshold: f64,
    #[case] expected_outliers: Vec<String>,
  ) -> Result<(), Report> {
    let (graph, names, times, branch_lengths) = setup_graph(IQD_TREE, &iqd_dates())?;
    let clock_model = ClockModel::for_testing(1.0, -2000.0);

    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let result = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, threshold, &NoopProgress)?;

    pretty_assert_abs_diff_eq!(2.0, result.iqd, epsilon = 1e-15);
    assert_eq!(expected_outliers, outlier_names(&names, &result.outliers));
    Ok(())
  }

  #[test]
  fn test_clock_filter_zero_interquartile_distance_flags_every_nonzero_deviation() -> Result<(), Report> {
    let dates = btreemap! {
      o!("A") => 2002.0,
      o!("B") => 2002.0,
      o!("C") => 2002.0,
      o!("D") => 2002.0,
      o!("E") => 2002.5,
    };
    let (graph, names, times, branch_lengths) = setup_graph("((A:1,B:1,C:1)X:1,(D:1,E:1)Y:1)root;", &dates)?;
    let clock_model = ClockModel::for_testing(1.0, -2000.0);

    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let result = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress)?;

    pretty_assert_abs_diff_eq!(0.0, result.iqd, epsilon = 1e-15);
    assert_eq!(vec![o!("E")], outlier_names(&names, &result.outliers));
    Ok(())
  }

  #[test]
  fn test_clock_filter_rejects_no_dated_leaves() -> Result<(), Report> {
    let (graph, times, branch_lengths) = setup_low_cardinality_graph(0)?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let result = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress);

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
    let (graph, times, branch_lengths) = setup_low_cardinality_graph(dated_leaf_count)?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let result = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress)?;

    assert!(result.iqd.is_finite());
    Ok(())
  }

  #[test]
  fn test_clock_filter_no_outliers_clean_data() -> Result<(), Report> {
    let (graph, _, times, branch_lengths) = setup_graph(helpers::FOUR_LEAF_TREE, &helpers::four_leaf_dates(2010.0))?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let result = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress)?;

    assert!(result.outliers.is_empty(), "No outliers expected for clean data");
    assert!(result.iqd.is_finite(), "IQD should be a finite number");
    assert!(result.iqd >= 0.0, "IQD should be non-negative");
    Ok(())
  }

  #[test]
  fn test_clock_filter_detects_outlier() -> Result<(), Report> {
    let (graph, names, times, branch_lengths) =
      setup_graph(helpers::FOUR_LEAF_TREE, &helpers::four_leaf_dates(1900.0))?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
    let result = clock_filter(&graph, &inputs, &clock_model, &branch_lengths, 3.0, &NoopProgress)?;

    assert!(result.iqd > 0.0, "IQD should be positive with varying dates");
    assert!(
      outlier_names(&names, &result.outliers).contains(&o!("A")),
      "Node A should be marked as outlier"
    );
    Ok(())
  }

  #[test]
  fn test_clock_filter_respects_threshold() -> Result<(), Report> {
    let (graph, _, times, branch_lengths) = setup_graph(helpers::FOUR_LEAF_TREE, &helpers::four_leaf_dates(1980.0))?;
    let clock_model = ClockModel::for_testing(0.01, -20.0);

    let inputs = ClockInputs::from_times(&graph, &times, &BTreeMap::new());
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

  mod helpers {
    use crate::o;
    use eyre::Report;
    use itertools::Itertools;
    use maplit::btreemap;
    use std::collections::{BTreeMap, BTreeSet};
    use treetime_graph::edge::GraphEdgeKey;
    use treetime_graph::graph::Graph;
    use treetime_graph::node::GraphNodeKey;
    use treetime_io::nwk::nwk_read_str;

    pub(super) const OUTLIER_TREE: &str = "(((A:0.1,B:0.2)AB:0.01,(C:0.15,D:0.25)CD:0.01)ABCD:0.01,((E:0.12,F:0.18)EF:0.01,(G:2.0,H:3.0)GH:0.01)EFGH:0.01)root;";

    pub(super) const IQD_TREE: &str = "((A:1,B:2,C:3,D:4,E:5)X:1,(F:6,G:7,H:8,I:9,J:100)Y:1)root;";

    pub(super) const FOUR_LEAF_TREE: &str = "((A:0.1,B:0.2)AB:0.1,(C:0.2,D:0.12)CD:0.05)root:0.01;";

    pub(super) type GraphSetup = (
      Graph,
      BTreeMap<GraphNodeKey, Option<String>>,
      BTreeMap<GraphNodeKey, Option<f64>>,
      BTreeMap<GraphEdgeKey, Option<f64>>,
    );

    pub(super) fn outlier_dates() -> BTreeMap<String, f64> {
      btreemap! {
        o!("A") => 2012.0,
        o!("B") => 2022.0,
        o!("C") => 2017.0,
        o!("D") => 2027.0,
        o!("E") => 2014.0,
        o!("F") => 2020.0,
        o!("G") => 2015.0,
        o!("H") => 2015.0,
      }
    }

    pub(super) fn outlier_tree_divergences() -> BTreeMap<String, f64> {
      btreemap! {
        o!("root") => 0.0,
        o!("ABCD") => 0.01,
        o!("EFGH") => 0.01,
        o!("AB") => 0.02,
        o!("CD") => 0.02,
        o!("EF") => 0.02,
        o!("GH") => 0.02,
        o!("A") => 0.12,
        o!("B") => 0.22,
        o!("C") => 0.17,
        o!("D") => 0.27,
        o!("E") => 0.14,
        o!("F") => 0.2,
        o!("G") => 2.02,
        o!("H") => 3.02,
      }
    }

    pub(super) fn iqd_dates() -> BTreeMap<String, f64> {
      btreemap! {
        o!("A") => 1996.0,
        o!("B") => 2002.0,
        o!("C") => 2004.0,
        o!("D") => 2005.5,
        o!("E") => 2007.0,
        o!("F") => 2008.5,
        o!("G") => 2010.0,
        o!("H") => 2012.0,
        o!("I") => 2016.5,
      }
    }

    pub(super) fn four_leaf_dates(date_a: f64) -> BTreeMap<String, f64> {
      btreemap! {
        o!("A") => date_a,
        o!("B") => 2020.0,
        o!("C") => 2015.0,
        o!("D") => 2012.0,
      }
    }

    pub(super) fn setup_graph(newick: &str, dates: &BTreeMap<String, f64>) -> Result<GraphSetup, Report> {
      let nwk_parsed = nwk_read_str(newick)?;
      let names = nwk_parsed.names();
      let graph = nwk_parsed.graph;
      let times = graph
        .get_leaves()
        .map(|leaf| {
          let name = names[&leaf.key()].clone().expect("leaves are named");
          (leaf.key(), dates.get(&name).copied())
        })
        .collect();
      Ok((graph, names, times, nwk_parsed.branch_lengths))
    }

    pub(super) fn setup_low_cardinality_graph(
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
      let graph = nwk_parsed.graph;
      let times = graph
        .get_leaves()
        .zip([2000.0, 2001.0, 2002.0])
        .take(dated_leaf_count)
        .map(|(leaf, date)| (leaf.key(), Some(date)))
        .collect();
      Ok((graph, times, nwk_parsed.branch_lengths))
    }

    pub(super) fn outlier_names(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      outliers: &BTreeSet<GraphNodeKey>,
    ) -> Vec<String> {
      outliers
        .iter()
        .map(|key| names[key].clone().expect("outliers are named"))
        .sorted()
        .collect()
    }

    pub(super) fn divergences_by_name(
      names: &BTreeMap<GraphNodeKey, Option<String>>,
      divergences: &BTreeMap<GraphNodeKey, f64>,
    ) -> BTreeMap<String, f64> {
      divergences
        .iter()
        .map(|(key, div)| (names[key].clone().expect("nodes are named"), *div))
        .collect()
    }
  }
}
