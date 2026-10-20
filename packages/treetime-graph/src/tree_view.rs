use crate::edge::GraphEdgeKey;
use crate::graph::Graph;
use crate::node::GraphNodeKey;
use crate::schedule::{GraphSchedule, TreeSchedule};
use eyre::Report;

pub struct TreeView<'g> {
  tree: TreeSchedule<'g>,
}

impl<'g> TreeView<'g> {
  pub fn new(graph: &'g Graph) -> Result<Self, Report> {
    Ok(Self {
      tree: GraphSchedule::new(graph)?.into_tree()?,
    })
  }

  pub fn root(&self) -> GraphNodeKey {
    self.tree.root()
  }

  pub fn preorder(&self) -> &[GraphNodeKey] {
    self.tree.preorder()
  }

  pub fn parent(&self, key: GraphNodeKey) -> Option<(GraphNodeKey, GraphEdgeKey)> {
    self.tree.parent(key)
  }

  pub fn children(&self, key: GraphNodeKey) -> &[(GraphNodeKey, GraphEdgeKey)] {
    self.tree.children(key)
  }
}
