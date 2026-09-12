use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NodeCommentProvider;
use treetime_utils::o;

/// Supplies the `date` Newick and Nexus comment for timetree output from the committed node times
/// held in the `times` value map, keyed by node id.
///
/// The date is formatted to two decimals to match v0's Newick output; see
/// `kb/decisions/timetree-nwk-date-two-decimal-precision.md`. Only dated nodes carry a comment.
pub struct DateCommentProvider<'a> {
  times: &'a BTreeMap<GraphNodeKey, f64>,
}

impl<'a> DateCommentProvider<'a> {
  pub fn new(times: &'a BTreeMap<GraphNodeKey, f64>) -> Self {
    Self { times }
  }
}

impl NodeCommentProvider for DateCommentProvider<'_> {
  fn node_comments(&self, key: GraphNodeKey) -> Result<BTreeMap<String, String>, Report> {
    Ok(match self.times.get(&key) {
      Some(time) => btreemap! { o!("date") => format!("{time:.2}") },
      None => BTreeMap::new(),
    })
  }
}
