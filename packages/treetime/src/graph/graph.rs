use crate::graph::collections::{Digraph as FastDigraph, Graph as FastGraph};
use crate::graph::core::{Edge, Traverse};
use eyre::Report;
use std::fmt::{Debug, Display};
use std::io::Write;
use std::sync::Arc;

pub struct Graph<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  g: FastDigraph<usize, N, E>,
  idx: usize,
}

impl<N, E> Graph<N, E>
where
  N: Clone + Debug + Display + Sync + Send,
  E: Clone + Debug + Display + Sync + Send,
{
  pub fn new() -> Self {
    Self {
      g: FastDigraph::<usize, N, E>::new(),
      idx: 0,
    }
  }

  pub fn add_node(&mut self, node_payload: N) -> usize {
    let idx = self.idx;
    self.g.add_node(idx, node_payload);
    self.idx += 1;
    idx
  }

  pub fn add_edge(&mut self, from_node_idx: usize, to_node_idx: usize, edge_payload: E) {
    self.g.add_edge(from_node_idx, to_node_idx, edge_payload);
  }

  pub fn iter_breadth_first_forward<F>(&mut self, explorer: F)
  where
    F: Fn(&Arc<Edge<usize, N, E>>) -> Traverse + Sync + Send + Copy,
  {
    self.g.par_breadth_first(0, explorer);
  }

  pub fn print<W: Write>(&self, mut writer: W) -> Result<(), Report> {
    self.g.print_graph(&mut writer)?;
    Ok(())
  }
}
