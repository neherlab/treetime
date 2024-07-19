#[cfg(test)]
pub mod tests {
  #![allow(clippy::many_single_char_names)]
  use crate::graph::graph::tests::{TestEdge, TestNode};
  use crate::graph::graph::Graph;
  use crate::o;
  use eyre::Report;

  #[rustfmt::skip]
  pub fn get_example_tree() -> Result<Graph<TestNode, TestEdge>, Report> {
    let mut graph = Graph::<TestNode, TestEdge>::new();

    let r1 = graph.add_node(TestNode(o!("r1")));
    let  a = graph.add_node(TestNode(o!("a")));
    let  b = graph.add_node(TestNode(o!("b")));
    let  c = graph.add_node(TestNode(o!("c")));
    let  d = graph.add_node(TestNode(o!("d")));
    let  e = graph.add_node(TestNode(o!("e")));
    let  f = graph.add_node(TestNode(o!("f")));
    let  g = graph.add_node(TestNode(o!("g")));
    let  h = graph.add_node(TestNode(o!("h")));
    let  i = graph.add_node(TestNode(o!("i")));


    graph.add_edge(r1, b, TestEdge(0.0))?;
    graph.add_edge(r1, a, TestEdge(0.1))?;
    graph.add_edge(a, c,  TestEdge(0.2))?;
    graph.add_edge(a, f,  TestEdge(0.3))?;
    graph.add_edge(b, d,  TestEdge(0.4))?;
    graph.add_edge(c, g,  TestEdge(0.5))?;
    graph.add_edge(d, e,  TestEdge(0.6))?;
    graph.add_edge(e, h,  TestEdge(0.7))?;
    graph.add_edge(f, i,  TestEdge(0.8))?;

    graph.build()?;

    Ok(graph)
  }

  #[rustfmt::skip]
  pub fn get_example_graph() -> Result<Graph<TestNode, TestEdge>, Report> {
    let mut graph = Graph::<TestNode, TestEdge>::new();

    let r1 = graph.add_node(TestNode(o!("r1")));
    let  a = graph.add_node(TestNode(o!("a")));
    let  b = graph.add_node(TestNode(o!("b")));
    let  c = graph.add_node(TestNode(o!("c")));
    let  d = graph.add_node(TestNode(o!("d")));
    let  e = graph.add_node(TestNode(o!("e")));
    let  f = graph.add_node(TestNode(o!("f")));
    let  g = graph.add_node(TestNode(o!("g")));
    let  h = graph.add_node(TestNode(o!("h")));
    let  i = graph.add_node(TestNode(o!("i")));


    graph.add_edge(r1, b, TestEdge(0.9))?;
    graph.add_edge(r1, a, TestEdge(1.0))?;
    graph.add_edge(a, c, TestEdge(1.1))?;
    graph.add_edge(a, d, TestEdge(1.2))?;
    graph.add_edge(a, f, TestEdge(1.3))?;
    graph.add_edge(b, d, TestEdge(1.4))?;
    graph.add_edge(c, g, TestEdge(1.5))?;
    graph.add_edge(d, e, TestEdge(1.6))?;
    graph.add_edge(d, g, TestEdge(1.7))?;
    graph.add_edge(e, h, TestEdge(1.8))?;
    graph.add_edge(f, i, TestEdge(1.9))?;

    graph.build()?;

    Ok(graph)
  }

  #[rustfmt::skip]
  pub fn get_example_graph_with_multiple_roots() -> Result<Graph<TestNode, TestEdge>, Report>  {
    let mut graph = Graph::<TestNode, TestEdge>::new();

    let r1 = graph.add_node(TestNode(o!("r1")));
    let r2 = graph.add_node(TestNode(o!("r2")));
    let  a = graph.add_node(TestNode(o!("a")));
    let  b = graph.add_node(TestNode(o!("b")));
    let  c = graph.add_node(TestNode(o!("c")));
    let  d = graph.add_node(TestNode(o!("d")));
    let  e = graph.add_node(TestNode(o!("e")));
    let  f = graph.add_node(TestNode(o!("f")));
    let  g = graph.add_node(TestNode(o!("g")));
    let  h = graph.add_node(TestNode(o!("h")));
    let  i = graph.add_node(TestNode(o!("i")));
    let  j = graph.add_node(TestNode(o!("j")));
    let  k = graph.add_node(TestNode(o!("k")));
    let  l = graph.add_node(TestNode(o!("l")));
    let  m = graph.add_node(TestNode(o!("m")));
    let  n = graph.add_node(TestNode(o!("n")));
    let  o = graph.add_node(TestNode(o!("o")));

    graph.add_edge(r1, a, TestEdge(2.0))?;
    graph.add_edge(r2, b, TestEdge(2.1))?;
    graph.add_edge(a, c, TestEdge(2.2))?;
    graph.add_edge(a, f, TestEdge(2.3))?;
    graph.add_edge(a, d, TestEdge(2.4))?;
    graph.add_edge(b, d, TestEdge(2.5))?;
    graph.add_edge(d, g, TestEdge(2.6))?;
    graph.add_edge(b, e, TestEdge(2.7))?;
    graph.add_edge(e, h, TestEdge(2.8))?;
    graph.add_edge(c, g, TestEdge(2.9))?;
    graph.add_edge(f, i, TestEdge(3.0))?;
    graph.add_edge(f, j, TestEdge(3.1))?;
    graph.add_edge(d, k, TestEdge(3.2))?;
    graph.add_edge(e, l, TestEdge(3.3))?;
    graph.add_edge(b, m, TestEdge(3.4))?;
    graph.add_edge(r2, n, TestEdge(3.5))?;
    graph.add_edge(r2, o, TestEdge(3.6))?;
    graph.add_edge(e, o, TestEdge(3.7))?;
    graph.add_edge(m, k, TestEdge(3.8))?;

    graph.build()?;

    Ok(graph)
  }
}
