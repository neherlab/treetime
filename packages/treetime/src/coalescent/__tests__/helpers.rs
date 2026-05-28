use crate::clock::date_constraints::load_date_constraints;
use crate::partition::timetree::GraphTimetree;
use eyre::Report;
use maplit::btreemap;
use treetime_io::dates_csv::DateConstraint;
use treetime_io::nwk::nwk_read_str;
use treetime_utils::o;

pub const TREE_NWK: &str = "((leaf1:0.01,leaf2:0.01)internal1:0.01,leaf3:0.02)root:0.0;";

pub fn setup_graph() -> Result<GraphTimetree, Report> {
  let dates = btreemap! {
    o!("root") => Some(DateConstraint::exact(2000.0)),
    o!("internal1") => Some(DateConstraint::exact(2005.0)),
    o!("leaf1") => Some(DateConstraint::exact(2010.0)),
    o!("leaf2") => Some(DateConstraint::exact(2010.0)),
    o!("leaf3") => Some(DateConstraint::exact(2012.0)),
  };
  let graph: GraphTimetree = nwk_read_str(TREE_NWK)?;
  load_date_constraints(&dates, &graph)?;
  Ok(graph)
}
