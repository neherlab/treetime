use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use crate::io::file::create_file_or_stdout;
use crate::io::nwk::{EdgeToNwk, NodeToNwk, NwkWriteOptions, nwk_write_str};
use eyre::Report;
use itertools::Itertools;
use smart_default::SmartDefault;
use std::io::Write;
use std::path::Path;

#[derive(Clone, SmartDefault)]
pub struct NexWriteOptions {
  /// Format node weights keeping this many significant digits
  pub weight_significant_digits: Option<u8>,

  /// Format node weights keeping this many decimal digits
  pub weight_decimal_digits: Option<i8>,
}

pub fn nex_write_file<N, E, D>(
  filepath: impl AsRef<Path>,
  graph: &Graph<N, E, D>,
  options: &NexWriteOptions,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  let mut f = create_file_or_stdout(filepath)?;
  nex_write(&mut f, graph, options)?;
  writeln!(f)?;
  Ok(())
}

pub fn nex_write_string<N, E, D>(graph: &Graph<N, E, D>, options: &NexWriteOptions) -> Result<String, Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  let mut buf = Vec::new();
  nex_write(&mut buf, graph, options)?;
  Ok(String::from_utf8(buf)?)
}

pub fn nex_write<N, E, D>(w: &mut impl Write, graph: &Graph<N, E, D>, options: &NexWriteOptions) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  let n_leaves = graph.num_leaves();
  let leaf_names = graph
    .get_leaves()
    .iter()
    .filter_map(|n| {
      let n = n.read_arc().payload().read_arc();
      n.nwk_name().map(|s| s.as_ref().to_owned())
    })
    .join(" ");
  let nwk = nwk_write_str(
    graph,
    &NwkWriteOptions {
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
