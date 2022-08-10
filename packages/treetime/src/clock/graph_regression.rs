#![allow(non_upper_case_globals)]
#![allow(non_snake_case)]

use crate::clock::clock_graph::{ClockGraph, Node};
use crate::clock::graph_regression_policy::GraphNodeRegressionPolicy;
use crate::graph::graph::{GraphNodeBackward, GraphNodeForward, NodeEdgePair};
use approx::assert_ulps_ne;
use eyre::Report;
use itertools::Itertools;
use ndarray::{array, Array1, Array2, ArrayBase, Data, Ix1};
use ndarray_linalg::Inverse;
use ndarray_stats::CorrelationExt;
use num_traits::real::Real;
use num_traits::Float;

const tavgii: usize = 0;
const davgii: usize = 1;
const tsqii: usize = 2;
const dtavgii: usize = 3;
const dsqii: usize = 4;
const sii: usize = 5;

pub fn calculate_averages<P: GraphNodeRegressionPolicy>(graph: &mut ClockGraph) {
  graph.par_iter_breadth_first_backward(
    |GraphNodeBackward {
       is_root,
       is_leaf,
       key,
       payload: n,
       children,
     }| {
      if is_leaf {
        return;
      }
      let mut Q = Array1::<f64>::zeros(6);
      for NodeEdgePair { node: c, edge } in children {
        let c = c.read();
        let edge = edge.read();
        let tv = P::tip_value(&c);
        let bv = P::branch_value(&edge);
        let var = P::branch_variance(&c);
        Q += &propagate_averages(&c, tv, bv, var, false);
      }
      n.Q = Q;
    },
  );

  graph.par_iter_breadth_first_forward(
    |GraphNodeForward {
       key,
       payload: node,
       parents,
       ..
     }| {
      if node.is_root() {
        node.Qtot = node.Q.clone();
        return;
      }

      let mut O = Array1::<f64>::zeros(6);

      for NodeEdgePair { node: parent, edge } in &parents {
        let parent = parent.read();
        let edge = edge.read();
        let tv = P::tip_value(&parent);
        let bv = P::branch_value(&edge);
        let var = P::branch_variance(&parent);
        O += &propagate_averages(&parent, tv, bv, var, false);
      }

      for NodeEdgePair { node: parent, edge } in &parents {
        let parent = parent.read();
        let edge = edge.read();
        if !node.is_root() {
          let tv = P::tip_value(&parent);
          let bv = P::branch_value(&edge);
          let var = P::branch_variance(&parent);
          O += &propagate_averages(&parent, tv, bv, var, true);
        }
      }

      node.O = O;

      if !node.is_leaf() {
        let bv = {
          if parents.len() > 1 {
            unimplemented!("Multiple parent nodes are not supported yet");
          }

          let mut bv = 0.0;
          for NodeEdgePair { edge, .. } in &parents {
            let edge = edge.read();
            bv += P::branch_value(&edge);
          }
          bv
        };

        let tv = P::tip_value(node);
        let var = P::branch_variance(node);
        node.Qtot = &node.Q + &propagate_averages(node, tv, bv, var, true);
      }
    },
  );
}

/// Propagates means, variance, and covariances along a branch. Operates both towards the root and tips.
pub fn propagate_averages(n: &Node, tv: Option<f64>, bv: f64, var: f64, outgroup: bool) -> Array1<f64> {
  if n.is_leaf() && !outgroup {
    match tv {
      None => Array1::<f64>::zeros(6),
      Some(tv) => {
        if tv.is_infinite() || tv.is_nan() {
          Array1::<f64>::zeros(6)
        } else if var == 0.0 {
          Array1::<f64>::from_elem(6, f64::infinity())
        } else {
          array![
            tv / var,
            bv / var,
            tv.powf(2.0) / var,
            bv * tv / var,
            bv.powf(2.0) / var,
            1.0 / var
          ]
        }
      }
    }
  } else {
    let tmpQ = if outgroup { &n.O } else { &n.Q };
    let denom = 1.0 / (1.0 + var * tmpQ[sii]);
    array![
      tmpQ[tavgii] * denom,
      (tmpQ[davgii] + bv * tmpQ[sii]) * denom,
      tmpQ[tsqii] - var * tmpQ[tavgii].powf(2.0) * denom,
      tmpQ[dtavgii] + tmpQ[tavgii] * bv - var * tmpQ[tavgii] * (tmpQ[davgii] + bv * tmpQ[sii]) * denom,
      tmpQ[dsqii] + 2.0 * bv * tmpQ[davgii] + bv.powf(2.0) * tmpQ[sii]
        - var
          * (tmpQ[davgii].powf(2.0) + 2.0 * bv * tmpQ[davgii] * tmpQ[sii] + bv.powf(2.0) * tmpQ[sii].powf(2.0))
          * denom,
      tmpQ[sii] * denom,
    ]
  }
}

#[derive(Debug, Clone)]
pub struct BaseRegressionResult {
  pub slope: f64,
  pub intercept: f64,
  pub chisq: f64,
  pub hessian: Option<f64>,
  pub cov: Option<f64>,
}

impl Default for BaseRegressionResult {
  fn default() -> Self {
    Self {
      slope: 0.0,
      intercept: 0.0,
      chisq: 0.0,
      hessian: None,
      cov: None,
    }
  }
}

/// Calculates the regression coefficients for a given vector containing the averages of tip and
/// branch quantities
pub fn base_regression<S>(Q: &ArrayBase<S, Ix1>, slope: &Option<f64>) -> BaseRegressionResult
where
  S: Data<Elem = f64>,
{
  assert!(Q.iter().all(|x| x.is_finite()), "Non-finite values are invalid here");
  assert_ulps_ne!(Q[sii], 0.0);
  let V = Q[tsqii] - Q[tavgii].powf(2.0) / Q[sii];
  assert!(
    V > 0.0,
    "No variation in sampling dates! Please specify your clock rate explicitly. Variance was: {V}"
  );

  let (slope, only_intercept) = match slope {
    None => {
      let slope = (Q[dtavgii] - Q[tavgii] * Q[davgii] / Q[sii]) / V;
      (slope, false)
    }
    Some(slope) => (*slope, true),
  };

  let intercept = (Q[davgii] - Q[tavgii] * slope) / Q[sii];
  let chisq = if V > 0.0 {
    0.5 * (Q[dsqii] - Q[davgii].powf(2.0) / Q[sii] - (Q[dtavgii] - Q[davgii] * Q[tavgii] / Q[sii]).powf(2.0) / V)
  } else {
    0.5 * (Q[dsqii] - Q[davgii].powf(2.0) / Q[sii])
  };

  let (hessian, cov) = if !only_intercept {
    let hessian = array![[Q[tsqii], Q[tavgii]], [Q[tavgii], Q[sii]]];
    let cov = hessian.inv();
    (Some(hessian), Some(cov))
  } else {
    (None, None)
  };

  BaseRegressionResult {
    slope,
    intercept,
    chisq,
    hessian: None,
    cov: None,
  }
}

pub fn explained_variance<P>(graph: &mut ClockGraph) -> Result<f64, Report>
where
  P: GraphNodeRegressionPolicy,
{
  graph.par_iter_breadth_first_forward(
    |GraphNodeForward {
       is_root,
       is_leaf,
       key,
       payload: node,
       parents,
     }| {
      if node.is_leaf() {
        return;
      }

      if parents.len() > 1 {
        unimplemented!("Multiple parent nodes are not supported yet");
      }
      for NodeEdgePair { node: parent, edge } in parents {
        let parent = parent.read();
        let edge = edge.read();
        node.v += parent.v + P::branch_value(&edge);
      }
    },
  );

  let tips_vs_branches = graph
    .get_leaves()
    .iter()
    .flat_map(|leaf| {
      let leaf = leaf.read().payload();
      let leaf = leaf.read();
      let tip = P::tip_value(&leaf).expect("leaf node is expected to have tip value, but none found");
      let branch = leaf.v;
      [tip, branch]
    })
    .collect_vec();

  let tips_vs_branches = Array2::<f64>::from_shape_vec((graph.num_leaves(), 2), tips_vs_branches)?;

  let corr = tips_vs_branches.t().pearson_correlation()?;

  Ok(corr[[0, 1]])
}
