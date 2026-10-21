#[cfg(test)]
pub(crate) mod tests {
  use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
  use crate::coalescent::mean_field::LineageExpectation;
  use crate::coalescent::skyline::{SkylineParams, SkylineResult, optimize_skyline};
  use crate::progress::NoopProgress;
  use crate::test_utils::constraint_expectation;
  use crate::test_utils::dates_by_node;
  use eyre::Report;
  use maplit::btreemap;
  use std::collections::BTreeMap;
  use std::sync::Arc;
  use treetime_graph::graph::Graph;
  use treetime_io::dates_csv::DateConstraint;
  use treetime_io::nwk::nwk_read;
  use treetime_utils::o;

  const TREE_NWK: &str = "((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;";

  pub(crate) fn setup_graph() -> Result<(Graph, DateConstraints), Report> {
    let dates = btreemap! {
      o!("root") => Some(DateConstraint::exact(2000.0)),
      o!("internal1") => Some(DateConstraint::exact(2005.0)),
      o!("leaf1") => Some(DateConstraint::exact(2010.0)),
      o!("leaf2") => Some(DateConstraint::exact(2010.0)),
      o!("leaf3") => Some(DateConstraint::exact(2012.0)),
    };
    graph_with_dates(TREE_NWK, &dates)
  }

  pub(crate) fn expectation(graph: &Graph, constraints: &DateConstraints) -> Arc<LineageExpectation> {
    constraint_expectation(graph, constraints).unwrap()
  }

  pub(crate) fn constant_skyline(expectation: &LineageExpectation) -> Result<SkylineResult, Report> {
    optimize_skyline(
      &SkylineParams {
        n_points: 1,
        ..SkylineParams::default()
      },
      expectation,
      &NoopProgress,
    )
  }

  pub(crate) fn tc(result: &SkylineResult) -> f64 {
    result.tc_schedule.values()[0]
  }

  pub(crate) fn graph_with_dates(
    tree_nwk: &str,
    dates: &BTreeMap<String, Option<DateConstraint>>,
  ) -> Result<(Graph, DateConstraints), Report> {
    let nwk_parsed = nwk_read(tree_nwk.as_bytes())?;
    let names = nwk_parsed.names.clone();
    let graph = nwk_parsed.graph;
    let constraints = load_date_constraints(&dates_by_node(dates.clone(), &graph, &names), &graph, &NoopProgress)?;
    Ok((graph, constraints))
  }
}
