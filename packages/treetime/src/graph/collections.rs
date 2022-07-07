//=============================================================================
// TEMPLATE COLLECTIONS
//=============================================================================

//! This module offers the `Graph` trait, which allows user to create a graph
//! easily out of thier own desired container type or use one of the templates.
//!
use crate::graph::core::{
  connect, directed_breadth_traversal, directed_depth_traversal, disconnect, parallel_directed_breadth_traversal,
  parallel_undirected_breadth_traversal, undirected_breadth_traversal, undirected_depth_traversal, Edge, Empty, Node,
  Traverse,
};
use std::io::Write;
use std::{
  collections::HashMap,
  fmt::{Debug, Display},
  hash::Hash,
  sync::{Arc, Weak},
};

/// This trait can be used to easily create a graph from a desired container type.
/// User must specify graph direction, construction and nopde retrieval. Rest is
/// implemented by the trait.
pub trait Graph<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  /// Create a new graph.
  fn new() -> Self;

  /// Direction of the graph.
  fn directed() -> bool;

  /// Add a node to the graph.
  fn add_node(&mut self, key: K, data: N) -> bool;

  /// Get an atomic reference to a node. If node can't
  /// be found, returns None.
  fn get_node(&self, node: K) -> Option<Arc<Node<K, N, E>>>;

  /// Iterate nodes witha  closure.
  fn iter_nodes(&self, f: &mut dyn FnMut(Arc<Node<K, N, E>>));

  /// Count the nodes in the graph.
  fn node_count(&self) -> usize;

  // ========================================================================

  /// Add a new edge to the graph.
  fn add_edge(&mut self, source: K, target: K, data: E) -> bool {
    let s = self.get_node(source);
    let t = self.get_node(target);
    match s {
      Some(src) => match t {
        Some(trg) => {
          connect(&src, &trg, data);
          true
        }
        None => false,
      },
      None => false,
    }
  }

  /// Delete an edge from the graph.
  fn del_edge(&mut self, source: K, target: K) -> bool {
    let s = self.get_node(source);
    let t = self.get_node(target);
    match s {
      Some(src) => match t {
        Some(trg) => {
          disconnect(&src, &trg);
          true
        }
        None => false,
      },
      None => false,
    }
  }

  /// Get an edge if it exists.
  fn get_edge(&self, source: K, target: K) -> Option<Arc<Edge<K, N, E>>> {
    let s = self.get_node(source);
    let t = self.get_node(target);
    match s {
      Some(ss) => match t {
        Some(tt) => ss.find_outbound(&tt),
        None => None,
      },
      None => None,
    }
  }

  /// Count the number of edges in the graph.
  fn edge_count(&self) -> usize {
    let r: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    self.iter_nodes(&mut |n| {
      let o = r.load(std::sync::atomic::Ordering::Relaxed);
      r.store(o + n.outbound().len(), std::sync::atomic::Ordering::Relaxed);
    });
    r.load(std::sync::atomic::Ordering::Relaxed)
  }

  /// Approximate the size of the graph.
  fn size_of(&self) -> usize {
    (self.node_count() * std::mem::size_of::<Node<K, N, E>>())
      + (self.edge_count() * std::mem::size_of::<Edge<K, N, E>>())
  }

  /// Depth first traversal of the graph.
  fn depth_first<F>(&self, source: K, explorer: F) -> Option<Vec<Weak<Edge<K, N, E>>>>
  where
    F: Fn(&Arc<Edge<K, N, E>>) -> Traverse,
  {
    match self.get_node(source) {
      Some(s) => {
        if Self::directed() {
          directed_depth_traversal(&s, explorer)
        } else {
          undirected_depth_traversal(&s, explorer)
        }
      }
      None => None,
    }
  }

  /// Breadth first traversal of the graph.
  fn breadth_first<F>(&self, source: K, explorer: F) -> Option<Vec<Weak<Edge<K, N, E>>>>
  where
    F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
  {
    match self.get_node(source) {
      Some(s) => {
        if Self::directed() {
          directed_breadth_traversal(&s, explorer)
        } else {
          undirected_breadth_traversal(&s, explorer)
        }
      }
      None => None,
    }
  }

  /// Parallel breadth first traversal of the graph.
  fn par_breadth_first<F>(&self, source: K, explorer: F) -> Option<Vec<Weak<Edge<K, N, E>>>>
  where
    F: Fn(&Arc<Edge<K, N, E>>) -> Traverse + Sync + Send + Copy,
  {
    match self.get_node(source) {
      Some(s) => {
        if Self::directed() {
          parallel_directed_breadth_traversal(&s, explorer)
        } else {
          parallel_undirected_breadth_traversal(&s, explorer)
        }
      }
      None => None,
    }
  }

  /// Print graph nodes.
  fn print_nodes<W: Write>(&self, mut writer: W) {
    self.iter_nodes(&mut |node| writeln!(writer, "	{}", node).unwrap());
  }

  /// Print graph edges.
  fn print_edges<W: Write>(&self, mut writer: W) {
    let sign = if Self::directed() { "->" } else { "--" };
    self.iter_nodes(&mut |node| {
      for edge in node.outbound().iter() {
        writeln!(
          writer,
          "	{} {} {} [label = \"{}\"]",
          edge.source().key(),
          sign,
          edge.target().key(),
          edge.load()
        )
        .unwrap();
      }
    });
  }

  /// Print graph in .dot format.
  fn print_graph<W: Write>(&self, mut writer: W) -> std::io::Result<()> {
    let name = if Self::directed() { "digraph" } else { "graph" };
    writeln!(writer, "{} {{", name)?;
    self.print_nodes(&mut writer);
    self.print_edges(&mut writer);
    writeln!(writer, "}}")?;
    Ok(())
  }
}

/// Undirected graph with arbitrary edge values. Underlying container type is
/// a `HashMap` which gives us fast lookup by key-value.
pub struct Ungraph<K, N = Empty, E = Empty>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  nodes: HashMap<K, Arc<Node<K, N, E>>>,
}

impl<K, N, E> Graph<K, N, E> for Ungraph<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn new() -> Self {
    Self { nodes: HashMap::new() }
  }

  fn directed() -> bool {
    false
  }

  fn add_node(&mut self, key: K, data: N) -> bool {
    if self.nodes.contains_key(&key) {
      false
    } else {
      let node = Arc::new(Node::new(key.clone(), data));
      self.nodes.insert(key, node);
      true
    }
  }

  fn get_node(&self, node: K) -> Option<Arc<Node<K, N, E>>> {
    match self.nodes.get(&node) {
      Some(n) => Some(n.clone()),
      None => None,
    }
  }

  fn iter_nodes(&self, f: &mut dyn FnMut(Arc<Node<K, N, E>>)) {
    for (_, node) in self.nodes.iter() {
      f(node.clone());
    }
  }

  fn node_count(&self) -> usize {
    self.nodes.len()
  }
}

/// Directed graph with arbitrary edge values. Underlying container type is
/// a `HashMap` which gives us fast lookup by key-value.
pub struct Digraph<K, N = Empty, E = Empty>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  nodes: HashMap<K, Arc<Node<K, N, E>>>,
}

impl<K, N, E> Graph<K, N, E> for Digraph<K, N, E>
where
  K: Hash + Eq + Clone + Debug + Display + Sync + Send,
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  fn new() -> Self {
    Self { nodes: HashMap::new() }
  }

  fn directed() -> bool {
    true
  }

  fn add_node(&mut self, key: K, data: N) -> bool {
    if self.nodes.contains_key(&key) {
      false
    } else {
      let node = Arc::new(Node::new(key.clone(), data));
      self.nodes.insert(key, node);
      true
    }
  }

  fn get_node(&self, node: K) -> Option<Arc<Node<K, N, E>>> {
    self.nodes.get(&node).cloned()
  }

  fn iter_nodes(&self, f: &mut dyn FnMut(Arc<Node<K, N, E>>)) {
    for node in self.nodes.values() {
      f(Arc::clone(node));
    }
  }

  fn node_count(&self) -> usize {
    self.nodes.len()
  }
}
