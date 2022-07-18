use crate::graph::breadth_first::{
  directed_breadth_first_traversal_backward, directed_breadth_first_traversal_forward,
};
use crate::graph::edge::Edge;
use crate::graph::node::Node;
use eyre::Report;
use itertools::Itertools;
use parking_lot::RwLock;
use petgraph::visit::Walker;
use std::borrow::{Borrow, BorrowMut};
use std::fmt::{Debug, Display};
use std::io::Write;
use std::iter::Map;
use std::ops::DerefMut;
use std::slice::Iter;
use std::sync::Arc;

pub type SafeNode<N, E> = Arc<RwLock<Node<N, E>>>;

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

/// Represents graph node during safe traversal
pub struct GraphNodeSafe<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  pub is_root: bool,
  pub is_leaf: bool,
  pub key: usize,
  pub payload: Arc<RwLock<N>>,
  pub children: Vec<NodeEdgePair<N, E>>,
  pub parents: Vec<NodeEdgePair<N, E>>,
}

impl<N, E> GraphNodeSafe<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send + Weighted,
{
  pub fn from_node(node: &Arc<RwLock<Node<N, E>>>) -> Self {
    let node = node.read();
    let is_leaf = node.is_leaf();
    let is_root = node.is_root();
    let key = node.key();

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

    let children_edges = node.outbound();
    let children = children_edges
      .iter()
      .map(|child_edge| {
        let edge: Arc<RwLock<E>> = child_edge.payload();
        let node: Arc<RwLock<N>> = child_edge.target().read().payload();
        NodeEdgePair { node, edge }
      })
      .collect_vec();

    let payload = node.payload();

    Self {
      is_root,
      is_leaf,
      key,
      payload,
      children,
      parents,
    }
  }
}

pub trait Weighted {
  fn weight(&self) -> f64 {
    0.0
  }
}

#[derive(Debug)]
pub struct Graph<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send + Weighted,
{
  nodes: Vec<Arc<RwLock<Node<N, E>>>>,
  idx: usize,
  roots: Vec<usize>,
  leaves: Vec<usize>,
}

impl<N, E> Graph<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send + Weighted,
{
  pub const fn new() -> Self {
    Self {
      nodes: Vec::new(),
      idx: 0,
      roots: vec![],
      leaves: vec![],
    }
  }

  pub fn get_node(&self, index: usize) -> Option<Arc<RwLock<Node<N, E>>>> {
    self.nodes.get(index).map(Arc::clone)
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn for_each<T, F>(&self, f: &mut dyn FnMut(GraphNodeSafe<N, E>)) {
    self.nodes.iter().for_each(|node| f(GraphNodeSafe::from_node(&node)));
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn map<T, F>(&self, mut f: F) -> Vec<T>
  where
    F: FnMut(GraphNodeSafe<N, E>) -> T,
  {
    self
      .nodes
      .iter()
      .map(|node| f(GraphNodeSafe::from_node(&node)))
      .collect_vec()
  }

  /// Iterates nodes synchronously and in unspecified order
  pub fn filter_map<T, F>(&self, mut f: F) -> Vec<T>
  where
    F: FnMut(GraphNodeSafe<N, E>) -> Option<T>,
  {
    self
      .nodes
      .iter()
      .filter_map(|node| {
        let node = node.read();
        let is_leaf = node.is_leaf();
        let is_root = node.is_root();
        let key = node.key();

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

        let children_edges = node.outbound();
        let children = children_edges
          .iter()
          .map(|child_edge| {
            let edge: Arc<RwLock<E>> = child_edge.payload();
            let node: Arc<RwLock<N>> = child_edge.target().read().payload();
            NodeEdgePair { node, edge }
          })
          .collect_vec();

        let payload = node.payload();

        f(GraphNodeSafe {
          is_root,
          is_leaf,
          key,
          payload,
          children,
          parents,
        })
      })
      .collect_vec()
  }

  pub fn node_count(&self) -> usize {
    self.nodes.len()
  }

  pub fn get_roots(&self) -> Vec<Arc<RwLock<Node<N, E>>>> {
    self.roots.iter().filter_map(|idx| self.get_node(*idx)).collect_vec()
  }

  pub fn add_node(&mut self, node_payload: N) -> usize {
    let idx = self.idx;
    let node = Arc::new(RwLock::new(Node::new(idx, node_payload)));
    self.nodes.insert(idx, node);
    self.idx += 1;
    idx
  }

  /// Add a new edge to the graph.
  pub fn add_edge(&mut self, src_idx: usize, dst_idx: usize, edge_payload: E) {
    // TODO: handle errors properly

    let src_mtx = self
      .get_node(src_idx)
      .ok_or_else(|| format!("When adding a graph edge {src_idx}->{dst_idx}: Node {src_idx} not found."))
      .unwrap();

    let dst_mtx = self
      .get_node(dst_idx)
      .ok_or_else(|| format!("When adding a graph edge {src_idx}->{dst_idx}: Node {dst_idx} not found."))
      .unwrap();

    assert_ne!(
      src_idx, dst_idx,
      "When adding a graph edge {src_idx}->{dst_idx}: Attempted to connect node {src_idx} to itself."
    );

    let (src, dst) = (src_mtx.read(), dst_mtx.read());

    let connected = src
      .outbound()
      .iter()
      .any(|edge| edge.target().read().key() == dst.key());

    assert!(
      !connected,
      "When adding a graph edge {src_idx}->{dst_idx}: Nodes {src_idx} and {dst_idx} are already connected."
    );

    let new_edge = Arc::new(Edge::new(
      Arc::downgrade(&src_mtx),
      Arc::downgrade(&dst_mtx),
      edge_payload,
    ));

    src.outbound_mut().push(Arc::clone(&new_edge));
    dst.inbound_mut().push(Arc::downgrade(&new_edge));
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
      .enumerate()
      .filter_map(|(i, node)| node.read().is_root().then(|| i))
      .collect_vec();

    self.leaves = self
      .nodes
      .iter()
      .enumerate()
      .filter_map(|(i, node)| node.read().is_leaf().then(|| i))
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
      .iter_mut()
      .for_each(|node| node.write().mark_as_not_visited());
  }

  fn print_fake_edges<W: Write>(&self, mut writer: W, nodes: &[&SafeNode<N, E>]) {
    // Fake edges needed to align a set of nodes beautifully
    let fake_edges = nodes
      .iter()
      .map(|node| node.read().key())
      .chunks(2)
      .into_iter()
      .enumerate()
      .map(|(i, pair)| {
        let pair = pair.collect_vec();
        if pair.len() != 2 {
          return "".to_owned();
        }

        let left = pair[0];
        let right = pair[1];
        let weight = 1000 + i * 100;
        format!("      {left}-> {right} [style=invis, weight={weight}]")
      })
      .join("\n");

    writeln!(
      writer,
      "\n    // fake edges for alignment of nodes\n    {{\n      rank=same\n{fake_edges}\n    }}"
    )
    .unwrap();
  }

  fn print_node<W: Write>(&self, mut writer: W, node: &SafeNode<N, E>) {
    writeln!(
      writer,
      "    {} [label = \"{}\"]",
      node.read().key(),
      node.read().payload().read()
    )
    .unwrap();
  }

  /// Print graph nodes.
  fn print_nodes<W: Write>(&self, mut writer: W) {
    let roots = self.nodes.iter().filter(|node| node.read().is_root()).collect_vec();
    let leaves = self.nodes.iter().filter(|node| node.read().is_leaf()).collect_vec();
    let internal = self
      .nodes
      .iter()
      .filter(|node| {
        let node = node.read();
        !node.is_leaf() && !node.is_root()
      })
      .collect_vec();

    writeln!(writer, "\n  subgraph roots {{").unwrap();

    for node in &roots {
      self.print_node(&mut writer, node);
    }

    self.print_fake_edges(&mut writer, &roots);

    writeln!(writer, "  }}\n\n  subgraph internals {{").unwrap();

    for node in &internal {
      self.print_node(&mut writer, node);
    }

    writeln!(writer, "  }}\n\n  subgraph leaves {{").unwrap();

    for node in &leaves {
      self.print_node(&mut writer, node);
    }

    self.print_fake_edges(&mut writer, &leaves);

    writeln!(writer, "  }}").unwrap();
  }

  /// Print graph edges.
  fn print_edges<W: Write>(&self, mut writer: W) {
    self.nodes.iter().for_each(&mut |node: &Arc<RwLock<Node<N, E>>>| {
      for edge in node.read().outbound().iter() {
        let payload = edge.payload();
        let payload = payload.read();
        writeln!(
          writer,
          "  {} -> {} [xlabel = \"{}\", weight=\"{}\"]",
          edge.source().read().key(),
          edge.target().read().key(),
          payload,
          payload.weight()
        )
        .unwrap();
      }
    });
  }

  /// Print graph in .dot format.
  pub fn print_graph<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
    write!(
      writer,
      r#"
digraph Phylogeny {{
  graph [rankdir=LR, overlap=scale, splines=ortho, nodesep=1.0, ordering=out];
  edge  [overlap=scale];
  node  [shape=box];
"#
    )?;
    self.print_nodes(&mut writer);
    writeln!(writer)?;
    self.print_edges(&mut writer);
    writeln!(writer, "}}")?;
    Ok(())
  }
}
