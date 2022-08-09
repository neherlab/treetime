use crate::clock::clock_graph::{Node, NodeType};

/// Defines how to obtain node values required for regression
pub trait GraphNodeRegressionPolicy {
  fn tip_value(node: &Node) -> Option<f64>;

  fn branch_value(node: &Node) -> f64;

  fn branch_variance(node: &Node) -> f64;
}

/// Defines how to obtain node values required for regression during rerooting
pub struct GraphNodeRegressionPolicyReroot;

impl GraphNodeRegressionPolicy for GraphNodeRegressionPolicyReroot {
  #[inline]
  fn tip_value(node: &Node) -> Option<f64> {
    if node.is_leaf() && !node.bad_branch {
      node.raw_date_constraint
    } else {
      None
    }
  }

  #[inline]
  fn branch_value(node: &Node) -> f64 {
    node.mutation_length
  }

  #[inline]
  fn branch_variance(node: &Node) -> f64 {
    match node.node_type {
      NodeType::Leaf(_) => 1.0,
      _ => 0.0,
    }
  }
}
