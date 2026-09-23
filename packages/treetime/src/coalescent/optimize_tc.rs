use crate::coalescent::node_time::CoalescentNodeTimes;
use crate::coalescent::skyline::{SkylineParams, optimize_skyline};
use eyre::Report;
use treetime_graph::graph::Graph;
use treetime_primitives::LogLh;

pub fn optimize_tc(graph: &Graph, node_times: &CoalescentNodeTimes) -> Result<OptimizeTcResult, Report> {
  let result = optimize_skyline(
    graph,
    &SkylineParams {
      n_points: 1,
      ..SkylineParams::default()
    },
    node_times,
  )?;
  Ok(OptimizeTcResult {
    tc: result.tc_values[0],
    log_tc_variance: result.log_tc_variances[0],
    tc_lower_bound: result.tc_lower_bounds[0],
    tc_upper_bound: result.tc_upper_bounds[0],
    likelihood: result.log_likelihood,
  })
}

pub struct OptimizeTcResult {
  pub tc: f64,
  pub log_tc_variance: f64,
  pub tc_lower_bound: f64,
  pub tc_upper_bound: f64,
  pub likelihood: LogLh,
}
