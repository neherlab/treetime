use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::coalescent::skyline::{SkylineParams, SkylineResult, optimize_skyline};
use crate::timetree::timetree_state::TimetreeState;
use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNodeKey;
use treetime_io::dates_csv::DateConstraint;
use treetime_io::nwk::nwk_read_str;
use treetime_utils::o;

const TREE_NWK: &str = "((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;";

pub(crate) fn setup_graph() -> Result<(Graph, BTreeMap<GraphNodeKey, Option<String>>, DateConstraints), Report> {
  let dates = btreemap! {
    o!("root") => Some(DateConstraint::exact(2000.0)),
    o!("internal1") => Some(DateConstraint::exact(2005.0)),
    o!("leaf1") => Some(DateConstraint::exact(2010.0)),
    o!("leaf2") => Some(DateConstraint::exact(2010.0)),
    o!("leaf3") => Some(DateConstraint::exact(2012.0)),
  };
  let nwk_parsed = nwk_read_str(TREE_NWK)?;
  let names = nwk_parsed.names();
  let graph = nwk_parsed.graph;
  let graph: Graph = graph;
  let constraints = load_date_constraints(&dates, &graph, &names)?;
  Ok((graph, names, constraints))
}

pub(crate) fn coalescent_node_times(graph: &Graph, constraints: &DateConstraints) -> CoalescentNodeTimes {
  TimetreeState::seed_from_values(graph, constraints).coalescent_node_times()
}

pub(crate) fn constant_skyline(graph: &Graph, node_times: &CoalescentNodeTimes) -> Result<SkylineResult, Report> {
  optimize_skyline(
    graph,
    &SkylineParams {
      n_points: 1,
      ..SkylineParams::default()
    },
    node_times,
  )
}

pub(crate) fn tc(result: &SkylineResult) -> f64 {
  result.tc_schedule.values()[0]
}
