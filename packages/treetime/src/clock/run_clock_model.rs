use crate::clock::clock_graph::ClockGraph;
use crate::clock::graph_regression::{base_regression, explained_variance, BaseRegressionResult};
use crate::clock::graph_regression_policy::{GraphNodeRegressionPolicy, GraphNodeRegressionPolicyReroot};
use crate::make_error;
use eyre::Report;
use itertools::Itertools;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

pub struct RunClockModelParams {
  pub slope: f64,
}

pub fn run_clock_model<P>(graph: &mut ClockGraph, params: &RunClockModelParams) -> Result<Vec<RootToTipResult>, Report>
where
  P: GraphNodeRegressionPolicy,
{
  let clock_model = regression(graph, params)?;

  if !clock_model.regression.slope.is_finite() {
    return make_error!(
      "Clock rate estimation failed. If your data lacks temporal signal, please specify the rate explicitly!"
    );
  }

  Ok(get_root_to_tip_result(graph, &clock_model))
}

pub struct ClockModel {
  regression: BaseRegressionResult,
  r_val: f64,
}

impl ClockModel {
  #[inline]
  pub const fn clock_rate(&self) -> f64 {
    self.regression.slope
  }

  #[inline]
  pub fn numdate_from_dist2root(&self, dist2root: f64) -> f64 {
    (dist2root - self.regression.intercept) / self.clock_rate()
  }

  #[inline]
  pub fn clock_deviation(&self, numdate: f64, dist2root: f64) -> f64 {
    (self.numdate_from_dist2root(dist2root) - numdate) * self.clock_rate()
  }
}

/// Regress tip values against branch values
fn regression(graph: &mut ClockGraph, params: &RunClockModelParams) -> Result<ClockModel, Report> {
  let Q_root = {
    if graph.num_roots() > 1 {
      unimplemented!("Multiple roots are not supported yet");
    }
    let roots = graph.get_roots();
    let root = roots[0].read();
    &root.payload().read().Q.clone()
  };

  let regression = base_regression(Q_root, &Some(params.slope));
  let r_val = explained_variance::<GraphNodeRegressionPolicyReroot>(graph)?;

  Ok(ClockModel { regression, r_val })
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct RootToTipResult {
  name: String,
  date: f64,
  #[serde(rename = "root-to-tip distance")]
  root_to_tip_distance: f64,
  #[serde(rename = "clock-deviation")]
  clock_deviation: f64,
}

/// Get results of the root-to-tip clock inference.
fn get_root_to_tip_result(graph: &ClockGraph, clock_model: &ClockModel) -> Vec<RootToTipResult> {
  graph
    .get_nodes()
    .par_iter()
    .map(|node| {
      let node = node.read().payload();
      let node = node.read();
      let name = node.name.clone();

      let (date, clock_deviation) = match node.raw_date_constraint {
        // TODO: is the condition on `is_leaf` necessary here? Could internal nodes have date constraints sometimes?
        Some(raw_date_constraint) if node.is_leaf() => {
          let clock_deviation = clock_model.clock_deviation(raw_date_constraint, node.dist2root);
          (raw_date_constraint, clock_deviation)
        }
        _ => {
          // Dates of nodes that didn't have a specified date are inferred from the root-to-tip regression.
          let date = clock_model.numdate_from_dist2root(node.dist2root);
          (date, 0.0)
        }
      };

      RootToTipResult {
        name,
        date,
        root_to_tip_distance: node.dist2root,
        clock_deviation,
      }
    })
    .collect()
}
