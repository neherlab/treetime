use eyre::Report;
use maplit::btreemap;
use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;
use treetime_io::nwk::NodeCommentProvider;
use treetime_utils::o;

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
