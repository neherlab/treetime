use std::collections::BTreeMap;
use treetime_graph::node::GraphNodeKey;

pub(crate) type CoalescentNodeTimes = BTreeMap<GraphNodeKey, CoalescentNodeTime>;

#[derive(Clone, Copy, Debug, Default)]
pub struct CoalescentNodeTime {
  pub time: Option<f64>,
  pub time_dist_likely: Option<f64>,
  pub bad_branch: bool,
}
