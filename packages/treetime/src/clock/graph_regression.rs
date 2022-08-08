#![allow(non_upper_case_globals)]

use crate::clock::clock_graph::{ClockGraph, Node, NodeType};
use crate::clock::graph_regression_policy::GraphNodeRegressionPolicy;
use crate::clock::run_reroot::RerootParams;
use crate::graph::graph::{GraphNodeBackward, NodeEdgePair};
use approx::{assert_ulps_ne, UlpsEq};
use ndarray::{array, Array1, ArrayBase, Data, Ix1};
use ndarray_linalg::Inverse;
use num_traits::real::Real;
use num_traits::Float;
use petgraph::matrix_graph::Zero;

const tavgii: usize = 0;
const davgii: usize = 1;
const tsqii: usize = 2;
const dtavgii: usize = 3;
const dsqii: usize = 4;
const sii: usize = 5;

pub fn calculate_averages<P: GraphNodeRegressionPolicy>(graph: &mut ClockGraph, params: &RerootParams) {
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
        let tv = P::tip_value(n);
        let bv = P::branch_value(n);
        let var = P::branch_variance(n);
        Q += &propagate_averages(&c, tv, bv, var, false);
      }
      n.Q = Q;
    },
  );

  graph.par_iter_breadth_first_backward(
    |GraphNodeBackward {
       is_root,
       is_leaf,
       key,
       payload: n,
       children,
     }| {
      let O = Array1::<f64>::zeros(6);
      if is_root {
        n.Qtot = n.Q.clone();
        return;
      }

      for NodeEdgePair { node: c, edge } in children {
        let tip_value = P::tip_value(n);
        let branch_value = P::branch_value(n);
        let branch_variance = P::branch_variance(n);
        // O += propagate_averages(c, tv, bv, var)

        // if n.up!=self.tree.root:
        //     c = n.up
        //     tv = self.tip_value(c)
        //     bv = self.branch_value(c)
        //     var = self.branch_variance(c)
        //     O += self.propagate_averages(c, tv, bv, var, outgroup=True)
        // n.O = O

        // if not n.is_terminal():
        //     tv = self.tip_value(n)
        //     bv = self.branch_value(n)
        //     var = self.branch_variance(n)
        //     n.Qtot = n.Q + self.propagate_averages(n, tv, bv, var, outgroup=True)
      }
    },
  );
}

/// Propagates means, variance, and covariances along a branch. Operates both towards the root and tips.
pub fn propagate_averages(n: &Node, tv: Option<f64>, bv: f64, var: f64, outgroup: bool) -> Array1<f64> {
  match n.node_type {
    NodeType::Leaf(_) if !outgroup => match tv {
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
    },
    _ => {
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
