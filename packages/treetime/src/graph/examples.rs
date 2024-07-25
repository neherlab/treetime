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

    let r1 = graph.add_node(TestNode(Some(o!("r1"))));
    let  a = graph.add_node(TestNode(Some(o!("a"))));
    let  b = graph.add_node(TestNode(Some(o!("b"))));
    let  c = graph.add_node(TestNode(Some(o!("c"))));
    let  d = graph.add_node(TestNode(Some(o!("d"))));
    let  e = graph.add_node(TestNode(Some(o!("e"))));
    let  f = graph.add_node(TestNode(Some(o!("f"))));
    let  g = graph.add_node(TestNode(Some(o!("g"))));
    let  h = graph.add_node(TestNode(Some(o!("h"))));
    let  i = graph.add_node(TestNode(Some(o!("i"))));


    graph.add_edge(r1, b, TestEdge(Some(0.0)))?;
    graph.add_edge(r1, a, TestEdge(Some(0.1)))?;
    graph.add_edge(a, c,  TestEdge(Some(0.2)))?;
    graph.add_edge(a, f,  TestEdge(Some(0.3)))?;
    graph.add_edge(b, d,  TestEdge(Some(0.4)))?;
    graph.add_edge(c, g,  TestEdge(Some(0.5)))?;
    graph.add_edge(d, e,  TestEdge(Some(0.6)))?;
    graph.add_edge(e, h,  TestEdge(Some(0.7)))?;
    graph.add_edge(f, i,  TestEdge(Some(0.8)))?;

    graph.build()?;

    Ok(graph)
  }

  #[rustfmt::skip]
  pub fn get_example_graph() -> Result<Graph<TestNode, TestEdge>, Report> {
    let mut graph = Graph::<TestNode, TestEdge>::new();

    let r1 = graph.add_node(TestNode(Some(o!("r1"))));
    let  a = graph.add_node(TestNode(Some(o!("a"))));
    let  b = graph.add_node(TestNode(Some(o!("b"))));
    let  c = graph.add_node(TestNode(Some(o!("c"))));
    let  d = graph.add_node(TestNode(Some(o!("d"))));
    let  e = graph.add_node(TestNode(Some(o!("e"))));
    let  f = graph.add_node(TestNode(Some(o!("f"))));
    let  g = graph.add_node(TestNode(Some(o!("g"))));
    let  h = graph.add_node(TestNode(Some(o!("h"))));
    let  i = graph.add_node(TestNode(Some(o!("i"))));


    graph.add_edge(r1, b, TestEdge(Some(0.9)))?;
    graph.add_edge(r1, a, TestEdge(Some(1.0)))?;
    graph.add_edge(a, c, TestEdge(Some(1.1)))?;
    graph.add_edge(a, d, TestEdge(Some(1.2)))?;
    graph.add_edge(a, f, TestEdge(Some(1.3)))?;
    graph.add_edge(b, d, TestEdge(Some(1.4)))?;
    graph.add_edge(c, g, TestEdge(Some(1.5)))?;
    graph.add_edge(d, e, TestEdge(Some(1.6)))?;
    graph.add_edge(d, g, TestEdge(Some(1.7)))?;
    graph.add_edge(e, h, TestEdge(Some(1.8)))?;
    graph.add_edge(f, i, TestEdge(Some(1.9)))?;

    graph.build()?;

    Ok(graph)
  }

  #[rustfmt::skip]
  pub fn get_example_graph_with_multiple_roots() -> Result<Graph<TestNode, TestEdge>, Report>  {
    let mut graph = Graph::<TestNode, TestEdge>::new();

    let r1 = graph.add_node(TestNode(Some(o!("r1"))));
    let r2 = graph.add_node(TestNode(Some(o!("r2"))));
    let  a = graph.add_node(TestNode(Some(o!("a"))));
    let  b = graph.add_node(TestNode(Some(o!("b"))));
    let  c = graph.add_node(TestNode(Some(o!("c"))));
    let  d = graph.add_node(TestNode(Some(o!("d"))));
    let  e = graph.add_node(TestNode(Some(o!("e"))));
    let  f = graph.add_node(TestNode(Some(o!("f"))));
    let  g = graph.add_node(TestNode(Some(o!("g"))));
    let  h = graph.add_node(TestNode(Some(o!("h"))));
    let  i = graph.add_node(TestNode(Some(o!("i"))));
    let  j = graph.add_node(TestNode(Some(o!("j"))));
    let  k = graph.add_node(TestNode(Some(o!("k"))));
    let  l = graph.add_node(TestNode(Some(o!("l"))));
    let  m = graph.add_node(TestNode(Some(o!("m"))));
    let  n = graph.add_node(TestNode(Some(o!("n"))));
    let  o = graph.add_node(TestNode(Some(o!("o"))));

    graph.add_edge(r1, a, TestEdge(Some(2.0)))?;
    graph.add_edge(r2, b, TestEdge(Some(2.1)))?;
    graph.add_edge(a, c, TestEdge(Some(2.2)))?;
    graph.add_edge(a, f, TestEdge(Some(2.3)))?;
    graph.add_edge(a, d, TestEdge(Some(2.4)))?;
    graph.add_edge(b, d, TestEdge(Some(2.5)))?;
    graph.add_edge(d, g, TestEdge(Some(2.6)))?;
    graph.add_edge(b, e, TestEdge(Some(2.7)))?;
    graph.add_edge(e, h, TestEdge(Some(2.8)))?;
    graph.add_edge(c, g, TestEdge(Some(2.9)))?;
    graph.add_edge(f, i, TestEdge(Some(3.0)))?;
    graph.add_edge(f, j, TestEdge(Some(3.1)))?;
    graph.add_edge(d, k, TestEdge(Some(3.2)))?;
    graph.add_edge(e, l, TestEdge(Some(3.3)))?;
    graph.add_edge(b, m, TestEdge(Some(3.4)))?;
    graph.add_edge(r2, n, TestEdge(Some(3.5)))?;
    graph.add_edge(r2, o, TestEdge(Some(3.6)))?;
    graph.add_edge(e, o, TestEdge(Some(3.7)))?;
    graph.add_edge(m, k, TestEdge(Some(3.8)))?;

    graph.build()?;

    Ok(graph)
  }
}
