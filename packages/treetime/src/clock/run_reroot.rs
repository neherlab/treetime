use crate::clock::clock_graph::{ClockGraph, Node};
use crate::clock::graph_regression_policy::{GraphNodeRegressionPolicy, GraphNodeRegressionPolicyReroot};
use crate::graph::graph::GraphNodeForward;
use crate::timetree::timetree_args::RerootMode;
use num_traits::Float;

pub struct RerootParams {
  pub reroot: RerootMode,
  pub slope: f64,
  pub force_positive: bool,
  pub keep_node_order: bool,
}

/// Determines the best root and reroots the tree to this value.
///
/// Note that this can change the parent child relations of the tree and values associated with
/// branches rather than nodes (e.g. confidence) might need to be re-evaluated afterwards
pub fn run_reroot(graph: &mut ClockGraph, params: &RerootParams) {
  let best_root = find_best_root::<GraphNodeRegressionPolicyReroot>(graph, params);

  // if best_root is None:
  //     raise ValueError("Rerooting failed!")
  // best_node = best_root["node"]
  //
  // x = best_root["split"]
  // if x<1e-5:
  //     new_node = best_node
  // elif x>1.0-1e-5:
  //     new_node = best_node.up
  // else:
  //     # create new node in the branch and root the tree to it
  //     new_node = Phylo.BaseTree.Clade()
  //
  //     # insert the new node in the middle of the branch
  //     # by simple re-wiring the links on the both sides of the branch
  //     # and fix the branch lengths
  //     new_node.branch_length = best_node.branch_length*(1-x)
  //     new_node.up = best_node.up
  //     new_node.clades = [best_node]
  //     new_node.up.clades = [k if k!=best_node else new_node
  //                           for k in best_node.up.clades]
  //
  //     best_node.branch_length *= x
  //     best_node.up = new_node
  //
  // new_node.rtt_regression = best_root
  // self.tree.root_with_outgroup(new_node)
  //
  // if not keep_node_order:
  //     self.tree.ladderize()
  // for n in self.tree.get_nonterminals(order='postorder'):
  //     for c in n:
  //         c.up=n
  // return best_root
}

#[derive(Debug, Clone)]
pub struct BestRoot {
  pub chisq: f64,
  pub cov: f64,
  pub hessian: f64,
  pub node: f64,
  pub split: f64,
}

impl Default for BestRoot {
  fn default() -> Self {
    Self {
      chisq: f64::infinity(),
      cov: 0.0,
      hessian: 0.0,
      node: 0.0,
      split: 0.0,
    }
  }
}

/// Determines the position on the tree that minimizes the bilinear product of the inverse
/// covariance and the data vectors.
fn find_best_root<P: GraphNodeRegressionPolicy>(graph: &mut ClockGraph, params: &RerootParams) -> BestRoot {
  // calculate_averages::<GraphNodeRegressionPolicyReroot>(graph, params);

  let best_root = BestRoot::default();

  graph.par_iter_breadth_first_forward(
    |GraphNodeForward {
       is_root,
       is_leaf,
       key,
       payload: n,
       parents,
     }| {
      if is_root {
        return;
      }

      let tv = P::tip_value(n);
      let bv = P::branch_value(n);
      let var = P::branch_variance(n);

      let (x, chisq) = find_optimal_root_along_branch(n, tv, bv, var, params);

      unimplemented!();

      // if chisq<best_root["chisq"]:
      //     tmpQ = self.propagate_averages(n, tv, bv*x, var*x) \
      //          + self.propagate_averages(n, tv, bv*(1-x), var*(1-x), outgroup=True)
      //     reg = base_regression(tmpQ, slope=slope)
      //     if reg["slope"]>=0 or (force_positive==False):
      //         best_root = {"node":n, "split":x}
      //         best_root.update(reg)
      // tip_value = lambda x:np.mean(x.raw_date_constraint) if (x.is_terminal() and (x.bad_branch is False)) else None
      // branch_value = lambda x:x.mutation_length
      // branch_variance = lambda x:1.0 if x.is_terminal() else 0.0
      // if is_leaf && !ancestral_args.reconstruct_tip_states {
      //   return;
      // }
      //
      // if parents.len() > 1 {
      //   unimplemented!("Multiple parent nodes not handled yet");
      // }
      //
      // let NodeEdgePair { edge, node: parent } = &parents[0];
      // let parent = parent.read();
      //
      // // choose the value of the Cx(i), corresponding to the state of the
      // // parent node i. This is the state of the current node
      // node.seq_ii = choose2(&parent.seq_ii, &node.joint_Cx.t());
      // node.seq = choose1(&node.seq_ii, &alphabet.alphabet);
    },
  );

  best_root
}

fn find_optimal_root_along_branch(n: &Node, tv: Option<f64>, bv: f64, var: f64, params: &RerootParams) -> (f64, f64) {
  unimplemented!();

  // from scipy.optimize import minimize_scalar
  // def chisq(x):
  //     tmpQ = self.propagate_averages(n, tv, bv*x, var*x) \
  //          + self.propagate_averages(n, tv, bv*(1-x), var*(1-x), outgroup=True)
  //     return base_regression(tmpQ, slope=slope)['chisq']
  //
  // if n.bad_branch or (n!=self.tree.root and n.up.bad_branch):
  //     return np.nan, np.inf
  //
  //
  // chisq_prox = np.inf if n.is_terminal() else base_regression(n.Qtot, slope=slope)['chisq']
  // chisq_dist = np.inf if n==self.tree.root else base_regression(n.up.Qtot, slope=slope)['chisq']
  //
  // grid = np.linspace(0.001,0.999,6)
  // chisq_grid = np.array([chisq(x) for x in grid])
  // min_chisq = chisq_grid.min()
  // if chisq_prox<=min_chisq:
  //     return 0.0, chisq_prox
  // elif chisq_dist<=min_chisq:
  //     return 1.0, chisq_dist
  // else:
  //     ii = np.argmin(chisq_grid)
  //     bounds = (0 if ii==0 else grid[ii-1], 1.0 if ii==len(grid)-1 else grid[ii+1])
  //     sol = minimize_scalar(chisq, bounds=bounds, method="bounded")
  //     if sol["success"]:
  //         return sol['x'], sol['fun']
  //     else:
  //         return np.nan, np.inf
}
