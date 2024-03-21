use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use crate::io::nwk::{write_nwk_str, WriteNwkOptions};
use eyre::Report;
use itertools::Itertools;
use smart_default::SmartDefault;
use std::io::Write;

#[derive(Clone, SmartDefault)]
pub struct WriteNexOptions {
  /// Format node weights keeping this many significant digits
  pub weight_significant_digits: Option<u8>,

  /// Format node weights keeping this many decimal digits
  pub weight_decimal_digits: Option<i8>,
}

pub fn write_nex<N, E>(w: &mut impl Write, graph: &Graph<N, E>, options: &WriteNexOptions) -> Result<(), Report>
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
  let nwk = write_nwk_str(
    graph,
    &WriteNwkOptions {
      weight_significant_digits: options.weight_significant_digits,
      weight_decimal_digits: options.weight_decimal_digits,
    },
  )?;

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

pub fn to_nex_string<N, E>(graph: &Graph<N, E>, options: &WriteNexOptions) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
{
  let mut buf = Vec::new();
  write_nex(&mut buf, graph, options)?;
  Ok(String::from_utf8(buf)?)
}
