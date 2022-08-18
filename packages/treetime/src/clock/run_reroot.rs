use crate::clock::clock_graph::{ClockGraph, Edge, Node, NodeType};
use crate::clock::graph_regression::{base_regression, propagate_averages};
use crate::clock::graph_regression_policy::{GraphNodeRegressionPolicy, GraphNodeRegressionPolicyReroot};
use crate::clock::minimize_scalar::minimize_scalar_brent_bounded;
use crate::graph::breadth_first::GraphTraversalContinuation;
use crate::graph::graph::{GraphNodeForward, NodeEdgePayloadPair};
use crate::graph::ladderize::ladderize;
use crate::graph::node::GraphNodeKey;
use crate::graph::reroot::reroot;
use crate::timetree::timetree_args::RerootMode;
use crate::{make_error, make_internal_report};
use eyre::Report;
use ndarray::Array1;
use ndarray_stats::QuantileExt;
use num_traits::Float;
use parking_lot::Mutex;
use std::sync::Arc;

pub struct RerootParams {
  pub reroot: RerootMode,
  pub slope: Option<f64>,
  /// only accept positive evolutionary rate estimates when rerooting the tree
  pub force_positive: bool,
  pub keep_node_order: bool,
}

/// Determines the best root and reroots the tree to this value.
///
/// Note that this can change the parent child relations of the tree and values associated with
/// branches rather than nodes (e.g. confidence) might need to be re-evaluated afterwards
pub fn run_reroot(graph: &mut ClockGraph, params: &RerootParams) -> Result<(), Report> {
  let best_root = match params.reroot {
    RerootMode::LeastSquares => find_best_root_least_squares::<GraphNodeRegressionPolicyReroot>(graph, params),
    RerootMode::MinDev => unimplemented!("'min dev' rerooting"),
    RerootMode::Oldest => unimplemented!("'oldest' rerooting"),
    RerootMode::ClockFilter => unimplemented!("'clock filter' rerooting"),
    RerootMode::Mrca => unimplemented!("'MRCA' rerooting"),
  }?;

  if let Some(best_root_node_key) = best_root.node_key {
    let best_node = graph
      .get_node(best_root_node_key)
      .ok_or_else(|| make_internal_report!("No node with index {}", best_root_node_key))?;

    let x = best_root.split;

    let new_root = if x < 1e-5 {
      best_node
    } else if x > 1.0 - 1e-5 {
      let parents = graph.parents_of(&best_node.read());
      if parents.len() > 1 {
        unimplemented!("Multiple parent nodes are not supported yet");
      }
      Arc::clone(&parents[0].0)
    } else {
      unimplemented!();
    };

    let old_root = {
      let roots = graph.get_roots();
      if roots.len() > 1 {
        unimplemented!("Multiple roots are not supported yet");
      }
      Arc::clone(&roots[0])
    };

    reroot(graph, &old_root, &new_root)?;

    // Recalculate some bookkeeping
    graph.build()?;

    if !params.keep_node_order {
      ladderize(graph);
    }

    Ok(())
  } else {
    make_error!("Rerooting failed: unable to find the best root node")
  }
}

#[derive(Debug, Clone)]
pub struct BestRoot {
  pub chisq: f64,
  pub node_key: Option<GraphNodeKey>,
  pub split: f64,
  pub hessian: Option<f64>,
  pub cov: Option<f64>,
}

impl Default for BestRoot {
  fn default() -> Self {
    Self {
      chisq: f64::infinity(),
      node_key: None,
      split: 0.0,
      cov: None,
      hessian: None,
    }
  }
}

/// Determine the node that, when the tree is rooted on this node, results
/// in the best regression of temporal constraints and root to tip distances.
///
/// Determines the position on the tree that minimizes the bilinear product of the inverse
/// covariance and the data vectors.
fn find_best_root_least_squares<P: GraphNodeRegressionPolicy>(
  graph: &mut ClockGraph,
  params: &RerootParams,
) -> Result<BestRoot, Report> {
  let best_root = {
    let best_root = Mutex::new(BestRoot::default());

    graph.par_iter_breadth_first_forward(
      |GraphNodeForward {
         is_root,
         is_leaf,
         key,
         payload: node,
         parents,
       }| {
        if is_root {
          return GraphTraversalContinuation::Continue;
        }

        let bv = {
          if parents.len() > 1 {
            unimplemented!("Multiple parent nodes are not supported yet");
          }

          let mut bv = 0.0;
          for (_, edge) in &parents {
            let edge = edge.read();
            bv += P::branch_value(&edge);
          }
          bv
        };

        let tv = P::tip_value(node);
        let var = P::branch_variance(node);

        let bad_branch = node.bad_branch || parents.iter().any(|(parent, _)| parent.read().bad_branch);

        let (x, chisq) = if bad_branch {
          (f64::nan(), f64::infinity())
        } else {
          find_optimal_root_along_branch(node, &parents, tv, bv, var, params).unwrap()
        };

        let mut best_root = best_root.lock();
        if chisq < best_root.chisq {
          let tmpQ = propagate_averages(node, tv, bv * x, var * x, false)
            + propagate_averages(node, tv, bv * (1.0 - x), var * (1.0 - x), true);
          let reg = base_regression(&tmpQ, &params.slope);

          if reg.slope >= 0.0 || !params.force_positive {
            *best_root = BestRoot {
              chisq: reg.chisq,
              cov: reg.cov,
              hessian: reg.hessian,
              node_key: Some(key),
              split: x,
            };
          }
        }

        GraphTraversalContinuation::Continue
      },
    );

    best_root.into_inner()
  };

  if best_root.node_key.is_none() {
    return make_error!("No valid root found!");
  }

  if let Some(hessian) = best_root.hessian {
    unimplemented!("calculate differentials with respect to x");
    // # calculate differentials with respect to x
    // deriv = []
    // n = best_root["node"]
    // tv = self.tip_value(n)
    // bv = self.branch_value(n)
    // var = self.branch_variance(n)
    // for dx in [-0.001, 0.001]:
    //     # y needs to be bounded away from 0 and 1 to avoid division by 0
    //     y = min(0.9999, max(0.0001, best_root["split"]+dx))
    //     tmpQ = self.propagate_averages(n, tv, bv*y, var*y) \
    //          + self.propagate_averages(n, tv, bv*(1-y), var*(1-y), outgroup=True)
    //     reg = base_regression(tmpQ, slope=slope)
    //     deriv.append([y,reg['chisq'], tmpQ[tavgii], tmpQ[davgii]])
    //
    // estimator_hessian = np.zeros((3,3))
    // estimator_hessian[:2,:2] = best_root['hessian']
    // estimator_hessian[2,2] = (deriv[0][1] + deriv[1][1] - 2.0*best_root['chisq'])/(deriv[0][0] - deriv[1][0])**2
    // # estimator_hessian[2,0] = (deriv[0][2] - deriv[1][2])/(deriv[0][0] - deriv[1][0])
    // # estimator_hessian[2,1] = (deriv[0][3] - deriv[1][3])/(deriv[0][0] - deriv[1][0])
    // estimator_hessian[0,2] = estimator_hessian[2,0]
    // estimator_hessian[1,2] = estimator_hessian[2,1]
    // best_root['hessian'] = estimator_hessian
    // best_root['cov'] = np.linalg.inv(estimator_hessian)
  }

  Ok(best_root)
}

fn find_optimal_root_along_branch(
  n: &Node,
  parents: &[NodeEdgePayloadPair<Node, Edge>],
  tv: Option<f64>,
  bv: f64,
  var: f64,
  params: &RerootParams,
) -> Result<(f64, f64), Report> {
  let chisq_prox = match n.node_type {
    NodeType::Leaf(_) => f64::infinity(),
    _ => base_regression(&n.Qtot, &params.slope).chisq,
  };

  let chisq_dist = if n.node_type == NodeType::Root {
    f64::infinity()
  } else {
    if parents.len() != 1 {
      unimplemented!("Multiple parent nodes are not supported yet");
    };
    let parent_node = &parents[0].0.read();
    base_regression(&parent_node.Qtot, &params.slope).chisq
  };

  let cost_function = |x: f64| {
    let tmpQ = propagate_averages(n, tv, bv * x, var * x, false)
      + propagate_averages(n, tv, bv * (1.0 - x), var * (1.0 - x), true);
    base_regression(&tmpQ, &params.slope).chisq
  };

  let grid = Array1::<f64>::linspace(0.001, 0.999, 6);
  let chisq_grid = Array1::from_iter(grid.map(|x| cost_function(*x)));
  let min_chisq = *chisq_grid.min()?;

  Ok(if chisq_prox <= min_chisq {
    (0.0, chisq_prox)
  } else if chisq_dist <= min_chisq {
    (1.0, chisq_dist)
  } else {
    let ii = chisq_grid.argmin()?;
    let bounds_from = if ii == 0 { 0.0 } else { grid[ii - 1] };
    let bounds_to = if ii == grid.len() - 1 { 1.0 } else { grid[ii + 1] };
    minimize_scalar_brent_bounded(cost_function, (bounds_from, bounds_to))
      .unwrap_or_else(|_| (f64::nan(), f64::infinity()))
  })
}