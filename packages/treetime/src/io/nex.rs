use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use crate::io::nwk::to_nwk_string;
use eyre::Report;
use itertools::Itertools;
use std::io::Write;

pub fn write_nex<N, E>(w: &mut impl Write, graph: &Graph<N, E>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let n_leaves = graph.num_leaves();
  let leaf_names = graph
    .get_leaves()
    .iter()
    .map(|n| {
      let n = n.read().payload();
      let n = n.read();
      n.name().to_owned()
    })
    .join(" ");
  let nwk = to_nwk_string(graph)?;

  writeln!(
    w,
    r#"#NEXUS
Begin Taxa;
  Dimensions NTax={n_leaves};
  TaxLabels {leaf_names};
End;
Begin Trees;
  Tree tree1={nwk};
End;
"#
  )?;

  Ok(())
}

pub fn to_nex_string<N, E>(graph: &Graph<N, E>) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let mut buf = Vec::new();
  write_nex(&mut buf, graph)?;
  Ok(String::from_utf8(buf)?)
}
