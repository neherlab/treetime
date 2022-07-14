use crate::graph::breadth_first::{
  directed_breadth_first_traversal_backward, directed_breadth_first_traversal_forward,
};
use crate::graph::edge::Edge;
use crate::graph::node::Node;
use eyre::Report;
use itertools::Itertools;
use parking_lot::RwLock;
use std::borrow::{Borrow, BorrowMut};
use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::io::Write;
use std::sync::Arc;

pub struct NodeEdgePair<N, E> {
  pub node: Arc<RwLock<N>>,
  pub edge: Arc<RwLock<E>>,
}

/// Represents graph node during forward traversal
pub struct GraphNodeForward<'n, N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: usize,
  pub payload: &'n mut N,
  pub parents: Vec<NodeEdgePair<N, E>>,
}

/// Represents graph node during backwards traversal
pub struct GraphNodeBackward<'n, N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: usize,
  pub payload: &'n mut N,
  pub children: Vec<NodeEdgePair<N, E>>,
}

pub struct Graph<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  nodes: HashMap<usize, Arc<RwLock<Node<N, E>>>>,
  idx: usize,
  roots: Vec<usize>,
  leaves: Vec<usize>,
}

impl<N, E> Graph<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  pub fn new() -> Self {
    Self {
      nodes: HashMap::new(),
      idx: 0,
      roots: vec![],
      leaves: vec![],
    }
  }

  pub fn get_node(&self, node: usize) -> Option<Arc<RwLock<Node<N, E>>>> {
    self.nodes.get(&node).cloned()
  }

  pub fn iter_nodes(&self, f: &mut dyn FnMut(Arc<RwLock<Node<N, E>>>)) {
    for node in self.nodes.values() {
      f(Arc::clone(node));
    }
  }

  pub fn node_count(&self) -> usize {
    self.nodes.len()
  }

  pub fn add_node(&mut self, node_payload: N) -> usize {
    let idx = self.idx;
    let node = Arc::new(RwLock::new(Node::new(idx, node_payload)));
    self.nodes.insert(idx, node);
    self.idx += 1;
    idx
  }

  /// Add a new edge to the graph.
  pub fn add_edge(&mut self, src_idx: usize, dst_idx: usize, edge_payload: E) -> bool {
    let src = self.get_node(src_idx);
    let dst = self.get_node(dst_idx);
    match (src, dst) {
      (Some(src_mtx), Some(dst_mtx)) => {
        let (src, dst) = (src_mtx.read(), dst_mtx.read());

        let connected = src
          .outbound()
          .iter()
          .any(|edge| edge.target().read().key() == dst.key());

        if !connected {
          let new_edge = Arc::new(Edge::new(
            Arc::downgrade(&src_mtx),
            Arc::downgrade(&dst_mtx),
            edge_payload,
          ));
          src.outbound_mut().push(Arc::clone(&new_edge));
          dst.inbound_mut().push(Arc::downgrade(&new_edge));
          true
        } else {
          false
        }
      }
      _ => false,
    }
  }

  /// Delete an edge from the graph.
  pub fn del_edge(&mut self, src: usize, dst: usize) -> bool {
    let src = self.get_node(src);
    let dst = self.get_node(dst);
    match (src, dst) {
      (Some(src_mtx), Some(dst_mtx)) => {
        let (src, dst) = (src_mtx.read(), dst_mtx.read());

        let mut idx: (usize, usize) = (0, 0);
        let mut flag = false;
        for (i, edge) in src.outbound().iter().enumerate() {
          if edge.target().read().key() == dst.key() {
            idx.0 = i;
            flag = true;
          }
        }
        for (i, edge) in dst.inbound().iter().enumerate() {
          if edge.upgrade().unwrap().source().read().key() == src.key() {
            idx.1 = i;
          }
        }
        if flag {
          src.outbound_mut().remove(idx.0);
          src.inbound_mut().remove(idx.1);
        }
        flag
      }
      _ => false,
    }
  }

  pub fn build(&mut self) -> Result<(), Report> {
    self.roots = self
      .nodes
      .iter()
      .filter_map(|(i, node)| node.read().is_root().then(|| *i))
      .collect_vec();

    self.leaves = self
      .nodes
      .iter()
      .filter_map(|(i, node)| node.read().is_leaf().then(|| *i))
      .collect_vec();

    Ok(())
  }

  pub fn par_iter_breadth_first_forward<F>(&mut self, explorer: F)
  where
    F: Fn(GraphNodeForward<N, E>) + Sync + Send,
  {
    let roots = self.roots.iter().filter_map(|idx| self.get_node(*idx)).collect_vec();

    directed_breadth_first_traversal_forward(roots.as_slice(), |node| {
      let is_leaf = node.is_leaf();
      let is_root = node.is_root();
      let key = node.key();

      let payload = node.payload();
      let mut payload = payload.write();
      let payload = payload.borrow_mut();

      let parent_edges = node.inbound();
      let parents = parent_edges
        .iter()
        .map(|parent_edge| {
          let parent_edge = parent_edge.upgrade().unwrap();
          let edge: Arc<RwLock<E>> = parent_edge.payload();
          let node: Arc<RwLock<N>> = parent_edge.source().read().payload();
          NodeEdgePair { node, edge }
        })
        .collect_vec();

      explorer(GraphNodeForward {
        is_root,
        is_leaf,
        key,
        payload,
        parents,
      });
    });

    self.reset_nodes();
  }

  pub fn par_iter_breadth_first_backward<F>(&mut self, explorer: F)
  where
    F: Fn(GraphNodeBackward<N, E>) + Sync + Send,
  {
    let leaves = self.leaves.iter().filter_map(|idx| self.get_node(*idx)).collect_vec();

    directed_breadth_first_traversal_backward(leaves.as_slice(), |node| {
      let is_leaf = node.is_leaf();
      let is_root = node.is_root();
      let key = node.key();

      let payload = node.payload();
      let mut payload = payload.write();
      let payload = payload.borrow_mut();

      let children_edges = node.outbound();
      let children = children_edges
        .iter()
        .map(|child_edge| {
          let edge: Arc<RwLock<E>> = child_edge.payload();
          let node: Arc<RwLock<N>> = child_edge.target().read().payload();
          NodeEdgePair { node, edge }
        })
        .collect_vec();

      explorer(GraphNodeBackward {
        is_root,
        is_leaf,
        key,
        payload,
        children,
      });
    });

    self.reset_nodes();
  }

  /// Returns graph into initial state after traversal
  pub fn reset_nodes(&mut self) {
    // Mark all nodes as not visited. As a part of traversal all nodes, one by one, marked as visited,
    // to ensure correctness of traversal. Here we reset the "is visited" markers this,
    // to allow for traversals again.
    self
      .nodes
      .values_mut()
      .for_each(|node| node.write().mark_as_not_visited());
  }

  /// Print graph nodes.
  fn print_nodes<W: Write>(&self, mut writer: W) {
    self.iter_nodes(&mut |node| writeln!(writer, "	{}", node.read()).unwrap());
  }

  /// Print graph edges.
  fn print_edges<W: Write>(&self, mut writer: W) {
    self.iter_nodes(&mut |node| {
      for edge in node.read().outbound().iter() {
        writeln!(
          writer,
          "	{} -> {} [label = \"{}\"]",
          edge.source().read().key(),
          edge.target().read().key(),
          edge.payload().read()
        )
        .unwrap();
      }
    });
  }

  /// Print graph in .dot format.
  pub fn print_graph<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
    writeln!(writer, "digraph {{")?;
    self.print_nodes(&mut writer);
    self.print_edges(&mut writer);
    writeln!(writer, "}}")?;
    Ok(())
  }
}
