use crate::graph::core::{
  connect, directed_breadth_first_traversal, disconnect, parallel_directed_breadth_first_traversal, Traverse,
};
use crate::graph::edge::Edge;
use crate::graph::node::Node;
use std::collections::HashMap;
use std::fmt::{Debug, Display};
use std::hash::Hash;
use std::io::Write;
use std::sync::{Arc, Weak};

pub struct Graph<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  nodes: HashMap<usize, Arc<Node<usize, N, E>>>,
  idx: usize,
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
    }
  }

  pub fn get_node(&self, node: usize) -> Option<Arc<Node<usize, N, E>>> {
    self.nodes.get(&node).cloned()
  }

  pub fn iter_nodes(&self, f: &mut dyn FnMut(Arc<Node<usize, N, E>>)) {
    for node in self.nodes.values() {
      f(Arc::clone(node));
    }
  }

  pub fn node_count(&self) -> usize {
    self.nodes.len()
  }

  pub fn add_node(&mut self, node_payload: N) -> usize {
    let idx = self.idx;
    let node = Arc::new(Node::new(idx, node_payload));
    self.nodes.insert(idx, node);
    self.idx += 1;
    idx
  }

  /// Add a new edge to the graph.
  pub fn add_edge(&mut self, src_idx: usize, dst_idx: usize, edge_payload: E) {
    let src = self.get_node(src_idx);
    let dst = self.get_node(dst_idx);
    match src {
      Some(src) => match dst {
        Some(dst) => {
          connect(&src, &dst, edge_payload);
        }
        None => {}
      },
      None => {}
    }
  }

  /// Delete an edge from the graph.
  pub fn del_edge(&mut self, source: usize, target: usize) -> bool {
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
  pub fn get_edge(&self, source: usize, target: usize) -> Option<Arc<Edge<usize, N, E>>> {
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
  pub fn edge_count(&self) -> usize {
    let r: std::sync::atomic::AtomicUsize = std::sync::atomic::AtomicUsize::new(0);
    self.iter_nodes(&mut |n| {
      let o = r.load(std::sync::atomic::Ordering::Relaxed);
      r.store(o + n.outbound().len(), std::sync::atomic::Ordering::Relaxed);
    });
    r.load(std::sync::atomic::Ordering::Relaxed)
  }

  /// Approximate the size of the graph.
  pub fn size_of(&self) -> usize {
    (self.node_count() * std::mem::size_of::<Node<usize, N, E>>())
      + (self.edge_count() * std::mem::size_of::<Edge<usize, N, E>>())
  }

  /// Breadth first traversal of the graph.
  pub fn breadth_first<F>(&self, source: usize, explorer: F) -> Option<Vec<Weak<Edge<usize, N, E>>>>
  where
    F: Fn(&Arc<Edge<usize, N, E>>) -> Traverse + Sync + Send + Copy,
  {
    match self.get_node(source) {
      Some(s) => directed_breadth_first_traversal(&s, explorer),
      None => None,
    }
  }

  pub fn par_iter_breadth_first_forward<F>(&mut self, explorer: F)
  where
    F: Fn(&Arc<Edge<usize, N, E>>) -> Traverse + Sync + Send + Copy,
  {
    self.par_breadth_first(0, explorer);
  }

  /// Parallel breadth first traversal of the graph.
  fn par_breadth_first<F>(&self, source: usize, explorer: F) -> Option<Vec<Weak<Edge<usize, N, E>>>>
  where
    F: Fn(&Arc<Edge<usize, N, E>>) -> Traverse + Sync + Send + Copy,
  {
    match self.get_node(source) {
      Some(s) => parallel_directed_breadth_first_traversal(&s, explorer),
      None => None,
    }
  }

  /// Print graph nodes.
  fn print_nodes<W: Write>(&self, mut writer: W) {
    self.iter_nodes(&mut |node| writeln!(writer, "	{}", node).unwrap());
  }

  /// Print graph edges.
  fn print_edges<W: Write>(&self, mut writer: W) {
    self.iter_nodes(&mut |node| {
      for edge in node.outbound().iter() {
        writeln!(
          writer,
          "	{} -> {} [label = \"{}\"]",
          edge.source().key(),
          edge.target().key(),
          edge.load()
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
