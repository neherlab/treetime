use crate::graph::edge::{GraphEdge, Weighted};
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named, NodeType, WithNwkComments};
use eyre::Report;
use std::fmt::{Display, Formatter};

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Node {
  pub name: String,
}

impl Named for Node {
  fn name(&self) -> &str {
    &self.name
  }

  fn set_name(&mut self, name: &str) {
    self.name = name.to_owned();
  }
}

impl GraphNode for Node {
  fn root(name: &str, weight: f64) -> Self {
    Self { name: name.to_owned() }
  }

  fn internal(name: &str, weight: f64) -> Self {
    Self { name: name.to_owned() }
  }

  fn leaf(name: &str) -> Self {
    Self { name: name.to_owned() }
  }

  fn set_node_type(&mut self, node_type: NodeType) {}
}

impl WithNwkComments for Node {}

impl Display for Node {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", self.name)
  }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Edge {
  pub name: String,
}

impl Weighted for Edge {
  fn weight(&self) -> f64 {
    1.0
  }
}

impl GraphEdge for Edge {
  fn new(weight: f64) -> Self {
    Self {
      name: format!("{weight}"),
    }
  }
}

impl Display for Edge {
  fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
    write!(f, "{}", &self.name)
  }
}

#[rustfmt::skip]
pub fn get_example_tree() -> Result<Graph<Node, Edge>, Report> {
  let mut graph = Graph::<Node, Edge>::new();

  let  r1 = graph.add_node(Node { name:  "r1".to_owned() });
  let  a = graph.add_node(Node { name:  "a".to_owned() });
  let  b = graph.add_node(Node { name:  "b".to_owned() });
  let  c = graph.add_node(Node { name:  "c".to_owned() });
  let  d = graph.add_node(Node { name:  "d".to_owned() });
  let  e = graph.add_node(Node { name:  "e".to_owned() });
  let  f = graph.add_node(Node { name:  "f".to_owned() });
  let  g = graph.add_node(Node { name:  "g".to_owned() });
  let  h = graph.add_node(Node { name:  "h".to_owned() });
  let  i = graph.add_node(Node { name:  "i".to_owned() });


  graph.add_edge(r1, b, Edge { name: "r1->b".to_owned() })?;
  graph.add_edge(r1, a, Edge { name: "r1->a".to_owned() })?;
  graph.add_edge(a,  c, Edge { name:  "a->c".to_owned() })?;
  graph.add_edge(a,  f, Edge { name:  "a->f".to_owned() })?;
  graph.add_edge(b,  d, Edge { name:  "b->d".to_owned() })?;
  graph.add_edge(c,  g, Edge { name:  "c->g".to_owned() })?;
  graph.add_edge(d,  e, Edge { name:  "d->e".to_owned() })?;
  graph.add_edge(e,  h, Edge { name:  "e->h".to_owned() })?;
  graph.add_edge(f,  i, Edge { name:  "f->i".to_owned() })?;

  graph.build()?;

  Ok(graph)
}

#[rustfmt::skip]
pub fn get_example_graph() -> Result<Graph<Node, Edge>, Report> {
  let mut graph = Graph::<Node, Edge>::new();

  let  r1 = graph.add_node(Node { name:  "r1".to_owned() });
  let  a = graph.add_node(Node { name:  "a".to_owned() });
  let  b = graph.add_node(Node { name:  "b".to_owned() });
  let  c = graph.add_node(Node { name:  "c".to_owned() });
  let  d = graph.add_node(Node { name:  "d".to_owned() });
  let  e = graph.add_node(Node { name:  "e".to_owned() });
  let  f = graph.add_node(Node { name:  "f".to_owned() });
  let  g = graph.add_node(Node { name:  "g".to_owned() });
  let  h = graph.add_node(Node { name:  "h".to_owned() });
  let  i = graph.add_node(Node { name:  "i".to_owned() });


  graph.add_edge(r1, b, Edge { name: "r1->b".to_owned() })?;
  graph.add_edge(r1, a, Edge { name: "r1->a".to_owned() })?;
  graph.add_edge(a,  c, Edge { name:  "a->c".to_owned() })?;
  graph.add_edge(a,  d, Edge { name:  "a->d".to_owned() })?;
  graph.add_edge(a,  f, Edge { name:  "a->f".to_owned() })?;
  graph.add_edge(b,  d, Edge { name:  "b->d".to_owned() })?;
  graph.add_edge(c,  g, Edge { name:  "c->g".to_owned() })?;
  graph.add_edge(d,  e, Edge { name:  "d->e".to_owned() })?;
  graph.add_edge(d,  g, Edge { name:  "d->g".to_owned() })?;
  graph.add_edge(e,  h, Edge { name:  "e->h".to_owned() })?;
  graph.add_edge(f,  i, Edge { name:  "f->i".to_owned() })?;

  graph.build()?;

  Ok(graph)
}

#[rustfmt::skip]
pub fn get_example_graph_with_multiple_roots() -> Result<Graph<Node, Edge>, Report>  {
  let mut graph = Graph::<Node, Edge>::new();

  let r1 = graph.add_node(Node { name: "r1".to_owned() });
  let r2 = graph.add_node(Node { name: "r2".to_owned() });
  let  a = graph.add_node(Node { name:  "a".to_owned() });
  let  b = graph.add_node(Node { name:  "b".to_owned() });
  let  c = graph.add_node(Node { name:  "c".to_owned() });
  let  d = graph.add_node(Node { name:  "d".to_owned() });
  let  e = graph.add_node(Node { name:  "e".to_owned() });
  let  f = graph.add_node(Node { name:  "f".to_owned() });
  let  g = graph.add_node(Node { name:  "g".to_owned() });
  let  h = graph.add_node(Node { name:  "h".to_owned() });
  let  i = graph.add_node(Node { name:  "i".to_owned() });
  let  j = graph.add_node(Node { name:  "j".to_owned() });
  let  k = graph.add_node(Node { name:  "k".to_owned() });
  let  l = graph.add_node(Node { name:  "l".to_owned() });
  let  m = graph.add_node(Node { name:  "m".to_owned() });
  let  n = graph.add_node(Node { name:  "n".to_owned() });
  let  o = graph.add_node(Node { name:  "o".to_owned() });

  graph.add_edge(r1, a, Edge { name: "r1->a".to_owned() })?;
  graph.add_edge(r2, b, Edge { name: "r2->b".to_owned() })?;
  graph.add_edge(a,  c, Edge { name:  "a->c".to_owned() })?;
  graph.add_edge(a,  f, Edge { name:  "a->f".to_owned() })?;
  graph.add_edge(a,  d, Edge { name:  "a->d".to_owned() })?;
  graph.add_edge(b,  d, Edge { name:  "b->d".to_owned() })?;
  graph.add_edge(d,  g, Edge { name:  "d->g".to_owned() })?;
  graph.add_edge(b,  e, Edge { name:  "b->e".to_owned() })?;
  graph.add_edge(e,  h, Edge { name:  "e->h".to_owned() })?;
  graph.add_edge(c,  g, Edge { name:  "c->g".to_owned() })?;
  graph.add_edge(f,  i, Edge { name:  "f->i".to_owned() })?;
  graph.add_edge(f,  j, Edge { name:  "f->j".to_owned() })?;
  graph.add_edge(d,  k, Edge { name:  "d->k".to_owned() })?;
  graph.add_edge(e,  l, Edge { name:  "e->l".to_owned() })?;
  graph.add_edge(b,  m, Edge { name:  "b->m".to_owned() })?;
  graph.add_edge(r2, n, Edge { name: "r2->n".to_owned() })?;
  graph.add_edge(r2, o, Edge { name: "r2->o".to_owned() })?;
  graph.add_edge(e,  o, Edge { name:  "e->o".to_owned() })?;
  graph.add_edge(m,  k, Edge { name:  "m->k".to_owned() })?;

  graph.build()?;

  Ok(graph)
}
