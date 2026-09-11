use crate::clock::date_constraints::{DateConstraints, load_date_constraints};
use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::partition::timetree::partition::GraphTimetree;
use crate::timetree::timetree_state::TimetreeState;
use eyre::Report;
use maplit::btreemap;
use treetime_graph::value_maps::node_names;
use treetime_io::dates_csv::DateConstraint;
use treetime_io::nwk::nwk_read_str;
use treetime_utils::o;

pub const TREE_NWK: &str = "((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;";

pub fn setup_graph() -> Result<(GraphTimetree, DateConstraints), Report> {
  let dates = btreemap! {
    o!("root") => Some(DateConstraint::exact(2000.0)),
    o!("internal1") => Some(DateConstraint::exact(2005.0)),
    o!("leaf1") => Some(DateConstraint::exact(2010.0)),
    o!("leaf2") => Some(DateConstraint::exact(2010.0)),
    o!("leaf3") => Some(DateConstraint::exact(2012.0)),
  };
  let graph: GraphTimetree = nwk_read_str(TREE_NWK)?.graph;
  let constraints = load_date_constraints(&dates, &graph, &node_names(&graph))?;
  Ok((graph, constraints))
}

/// Build the coalescent node-time value the collectors consume, from the date constraints
/// [`load_date_constraints`] returns, matching what the payload-reading seed produced right after the
/// constraints load: each node's committed time starts `None`, its distribution peak comes from the
/// constraint, and its bad-branch flag from the constraint pass.
pub fn coalescent_node_times(graph: &GraphTimetree, constraints: &DateConstraints) -> CoalescentNodeTimes {
  TimetreeState::seed_from_values(graph, constraints).coalescent_node_times()
}
