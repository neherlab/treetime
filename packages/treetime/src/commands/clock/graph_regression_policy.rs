use crate::commands::clock::clock_graph::{Edge, Node};

/// Defines how to obtain node values required for regression
pub trait GraphNodeRegressionPolicy {
  fn tip_value(node: &Node) -> Option<f64>;

  fn branch_value(edge: &Edge) -> f64;

  fn branch_variance(node: &Node) -> f64;
}

/// Defines how to obtain node values required for regression during rerooting
pub struct GraphNodeRegressionPolicyReroot;

impl GraphNodeRegressionPolicy for GraphNodeRegressionPolicyReroot {
  #[inline]
  fn tip_value(node: &Node) -> Option<f64> {
    // FIXME
    None

    // if node.is_leaf() && !node.bad_branch {
    //   node.raw_date_constraint
    // } else {
    //   None
    // }
  }

  #[inline]
  fn branch_value(edge: &Edge) -> f64 {
    // FIXME
    0.0
    // edge.weight
  }

  #[inline]
  fn branch_variance(_: &Node) -> f64 {
    // FIXME
    0.0

    // match node.node_type {
    //   NodeType::Leaf(_) => 1.0,
    //   _ => 0.0,
    // }
  }
}
