#[cfg(test)]
mod tests {
  use crate::reroot::params::{RerootMethod, RerootSpec};
  use eyre::Report;
  use helpers::{
    CLOCK_FILTER_IQD, DATES, OUTLIER_DATES, OUTLIER_TREE, TREE, UNDATED_TIP_DATES, UNDATED_TIP_TREE,
    date_free_root_child_lengths, root_child_lengths, run_clock,
  };
  use treetime_utils::pretty_assert_map_ulps_eq;

  #[test]
  fn test_pipeline_reroot_oldest_gives_the_same_root_with_and_without_the_outlier_filter() -> Result<(), Report> {
    let spec = RerootSpec::Method(RerootMethod::Oldest);

    let without_filter = run_clock(TREE, &DATES, &spec, 0.0)?;
    let with_filter = run_clock(TREE, &DATES, &spec, CLOCK_FILTER_IQD)?;

    pretty_assert_map_ulps_eq!(
      root_child_lengths(&without_filter)?,
      root_child_lengths(&with_filter)?,
      max_ulps = 4
    );
    Ok(())
  }

  #[test]
  fn test_pipeline_reroot_min_dev_roots_once_before_the_outlier_filter() -> Result<(), Report> {
    let spec = RerootSpec::Method(RerootMethod::MinDev);
    let expected = date_free_root_child_lengths(OUTLIER_TREE, &OUTLIER_DATES, &spec)?;

    let output = run_clock(OUTLIER_TREE, &OUTLIER_DATES, &spec, CLOCK_FILTER_IQD)?;

    assert!(
      !output.outliers.is_empty(),
      "the outlier filter must flag a tip in this case"
    );
    pretty_assert_map_ulps_eq!(expected, root_child_lengths(&output)?, max_ulps = 4);
    Ok(())
  }

  #[test]
  fn test_pipeline_reroot_min_dev_counts_an_undated_tip() -> Result<(), Report> {
    let spec = RerootSpec::Method(RerootMethod::MinDev);
    let expected = date_free_root_child_lengths(UNDATED_TIP_TREE, &UNDATED_TIP_DATES, &spec)?;

    let output = run_clock(UNDATED_TIP_TREE, &UNDATED_TIP_DATES, &spec, 0.0)?;

    pretty_assert_map_ulps_eq!(expected, root_child_lengths(&output)?, max_ulps = 4);
    Ok(())
  }

  mod helpers {
    use crate::cancel::NoopCancel;
    use crate::clock::pipeline::{self, ClockInput, ClockOutput, ClockParams};
    use crate::progress::NoopProgress;
    use crate::reroot::date_free::{LengthsAdvice, reroot_date_free};
    use crate::reroot::params::RerootSpec;
    use crate::test_utils::dates_by_node;
    use eyre::Report;
    use std::collections::BTreeMap;
    use treetime_graph::key_map::{EdgeMap, NodeMap};

    use treetime_graph::graph::Graph;

    use treetime_graph::node_names::NodeNames;
    use treetime_graph::schedule::GraphSchedule;
    use treetime_io::nwk::nwk_read;
    use treetime_primitives::date::DateConstraint;

    pub(super) const CLOCK_FILTER_IQD: f64 = 3.0;

    pub(super) const TREE: &str = "(((A:0.010,B:0.020)AB:0.010,(C:0.030,D:0.020)CD:0.015)X:0.020,\
                                   ((E:0.040,F:0.050)EF:0.010,(G:0.035,H:0.060)GH:0.020)Y:0.010)R;";

    pub(super) const DATES: [(&str, f64); 8] = [
      ("A", 2000.0),
      ("B", 2002.0),
      ("C", 2004.0),
      ("D", 2003.0),
      ("E", 2006.0),
      ("F", 2008.0),
      ("G", 2007.0),
      ("H", 2011.0),
    ];

    pub(super) const OUTLIER_TREE: &str = "((A:0.014,B:0.010,C:0.018,D:0.015,E:0.025,F:0.023,O:0.0495)X:0.010,\
                                           (P:0.006,Q:0.016,S:0.014,T:0.023,U:0.019,V:0.027)Y:0.010)R;";

    pub(super) const OUTLIER_DATES: [(&str, f64); 13] = [
      ("A", 2000.0),
      ("B", 2003.0),
      ("C", 2006.0),
      ("D", 2009.0),
      ("E", 2012.0),
      ("F", 2015.0),
      ("O", 2007.5),
      ("P", 2000.0),
      ("Q", 2003.0),
      ("S", 2006.0),
      ("T", 2009.0),
      ("U", 2012.0),
      ("V", 2015.0),
    ];

    pub(super) const UNDATED_TIP_TREE: &str =
      "((A:0.010,B:0.020)AB:0.010,(C:0.010,(D:0.020,U:0.300)DU:0.010)CDU:0.010)R;";

    pub(super) const UNDATED_TIP_DATES: [(&str, f64); 4] = [("A", 2000.0), ("B", 2002.0), ("C", 2001.0), ("D", 2003.0)];

    pub(super) fn run_clock(
      newick: &str,
      dates: &[(&str, f64)],
      spec: &RerootSpec,
      clock_filter: f64,
    ) -> Result<ClockOutput, Report> {
      let parsed = nwk_read(newick.as_bytes())?;
      let input = ClockInput {
        dates: node_dates(dates, &parsed.graph, &parsed.names),
        graph: parsed.graph,
        names: parsed.names,
        branch_lengths: parsed.branch_lengths,
      };
      let params = ClockParams {
        clock_filter,
        allow_negative_rate: true,
        reroot_spec: Some(spec.clone()),
      };
      Ok(pipeline::run(
        &params,
        input,
        &NoopCancel,
        &NoopProgress,
        &NoopProgress,
      )?)
    }

    pub(super) fn date_free_root_child_lengths(
      newick: &str,
      dates: &[(&str, f64)],
      spec: &RerootSpec,
    ) -> Result<BTreeMap<String, f64>, Report> {
      let parsed = nwk_read(newick.as_bytes())?;
      let mut names = parsed.names;
      let mut graph = parsed.graph;
      let mut branch_lengths = parsed.branch_lengths;
      let dated = node_dates(dates, &graph, &names);
      reroot_date_free(
        &mut graph,
        &mut branch_lengths,
        spec,
        &|key| dated.contains_key(key),
        LengthsAdvice::OptimizeTheTree,
        &mut names,
        &NoopProgress,
      )?;
      child_lengths(&graph, &names, &branch_lengths)
    }

    pub(super) fn root_child_lengths(output: &ClockOutput) -> Result<BTreeMap<String, f64>, Report> {
      child_lengths(&output.graph, &output.names, &output.branch_lengths)
    }

    fn node_dates(dates: &[(&str, f64)], graph: &Graph, names: &NodeNames) -> NodeMap<DateConstraint> {
      let dates = dates
        .iter()
        .map(|&(name, date)| (name.to_owned(), Some(DateConstraint::exact(date))));
      dates_by_node(dates, graph, names)
    }

    fn child_lengths(
      graph: &Graph,
      names: &NodeNames,
      branch_lengths: &EdgeMap<Option<f64>>,
    ) -> Result<BTreeMap<String, f64>, Report> {
      let schedule = GraphSchedule::new(graph)?;
      let root_key = graph.get_exactly_one_root()?.key();
      Ok(
        schedule
          .children(root_key)
          .iter()
          .map(|&(child_key, edge_key)| {
            let length = branch_lengths[&edge_key].expect("every edge has a length");
            (names[child_key].to_owned(), length)
          })
          .collect(),
      )
    }
  }
}
