use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;

pub(crate) type CoalescentNodeTimes = BTreeMap<GraphNodeKey, CoalescentNodeTime>;

#[derive(Clone, Copy, Debug, Default)]
pub struct CoalescentNodeTime {
  pub(crate) time: Option<f64>,
  pub(crate) time_dist_likely: Option<f64>,
  pub(crate) bad_branch: bool,
}
