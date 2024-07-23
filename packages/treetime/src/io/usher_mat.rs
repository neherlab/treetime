use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::io::file::create_file_or_stdout;
use crate::io::file::open_file_or_stdin;
use crate::io::nwk::{nwk_read_str, nwk_write_str, EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions};
use bytes::Buf;
use eyre::{Report, WrapErr};
use std::io::{Read, Write};
use std::path::Path;

pub use usher_mat_utils::{UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode};

pub fn usher_mat_pb_read_file<N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: UsherDataToGraphData + Sync + Send + Default,
  (): UsherToGraph<N, E, D>,
{
  let filepath = filepath.as_ref();
  usher_mat_pb_read(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading Usher MAT protobuf file '{filepath:#?}'"))
}

pub fn usher_mat_pb_read_bytes<N, E, D>(buf: impl Buf) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: UsherDataToGraphData + Sync + Send + Default,
  (): UsherToGraph<N, E, D>,
{
  let tree = usher_mat_utils::usher_mat_pb_read_bytes(buf).wrap_err("When reading Usher MAT protobuf bytes")?;
  usher_to_graph(&tree)
}

pub fn usher_mat_pb_read<N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: UsherDataToGraphData + Sync + Send + Default,
  (): UsherToGraph<N, E, D>,
{
  let tree = usher_mat_utils::usher_mat_pb_read(reader).wrap_err("When reading Usher MAT protobuf")?;
  usher_to_graph(&tree)
}

pub fn usher_mat_pb_write_file<N, E, D>(filepath: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  (): UsherFromGraph<N, E, D>,
{
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  usher_mat_pb_write(&mut f, graph)
    .wrap_err_with(|| format!("When reading Usher MAT protobuf file '{filepath:#?}'"))?;
  writeln!(f)?;
  Ok(())
}

pub fn usher_mat_pb_write_bytes<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  (): UsherFromGraph<N, E, D>,
{
  let tree = usher_from_graph(graph)?;
  let mut buf = Vec::new();
  usher_mat_utils::usher_mat_pb_write_bytes(&mut buf, &tree).wrap_err("When writing Usher MAT protobuf bytes")
}

pub fn usher_mat_pb_write<N, E, D>(writer: &mut impl Write, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  (): UsherFromGraph<N, E, D>,
{
  let tree = usher_from_graph(graph)?;
  usher_mat_utils::usher_mat_pb_write(writer, &tree).wrap_err("When writing Usher MAT protobuf")
}

/// Describes conversion from Usher tree global data when reading from Usher MAT protobuf
pub trait UsherDataToGraphData: Sized + Default {
  fn usher_data_to_graph_data(_: &UsherTree) -> Result<Self, Report> {
    Ok(Self::default())
  }
}

pub struct UsherTreeContext<'a> {
  pub node: &'a UsherTreeNode,
  pub mutations: &'a UsherMutationList,
  pub metadata: &'a UsherMetadata,
  pub tree: &'a UsherTree,
}

/// Describes conversion from Usher tree node data when reading from Usher MAT protobuf
pub trait UsherToGraph<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: UsherDataToGraphData + Sync + Send,
{
  fn usher_node_to_graph_components(context: &UsherTreeContext) -> Result<(N, E), Report>;
}

/// Convert Usher MAT protobuf to graph
pub fn usher_to_graph<N, E, D>(data: &UsherTree) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: UsherDataToGraphData + Sync + Send + Default,
  (): UsherToGraph<N, E, D>,
{
  println!("{}", &data.newick);

  let mut graph: Graph<N, E, D> = nwk_read_str(&data.newick)?;
  graph.set_data(D::usher_data_to_graph_data(data)?);

  // let mut i = 0;
  // graph.iter_depth_first_preorder_forward(|mut node| {
  //   let context = UsherTreeContext {
  //     node: &data.condensed_nodes[i],
  //     mutations: &data.node_mutations[i],
  //     metadata: &data.metadata[i],
  //     tree: data,
  //   };
  //
  //   // if let Some(name) = node.payload.name() {
  //   //   assert_eq!(context.node.node_name, name.as_ref());
  //   // }
  //
  //   let (graph_node, graph_edge) = <() as UsherToGraph<N, E, D>>::usher_node_to_graph_components(&context).unwrap();
  //
  //   *node.payload = graph_node;
  //
  //   // if let Some(parent) = node.parents.first() {
  //   //   let (_, edge) = parent;
  //   //   *edge.write_arc() = graph_edge;
  //   // }
  //
  //   i += 1;
  // });

  Ok(graph)
}

pub struct UsherGraphContext<'a, N, E, D>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  (): UsherFromGraph<N, E, D>,
{
  pub node: &'a N,
  pub edge: Option<&'a E>,
  pub graph: &'a Graph<N, E, D>,
}

/// Describes conversion to Usher tree node data when writing to Usher MAT protobuf
pub trait UsherFromGraph<N, E, D>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  (): UsherFromGraph<N, E, D>,
{
  fn usher_node_from_graph_components(
    context: &UsherGraphContext<N, E, D>,
  ) -> Result<(UsherTreeNode, UsherMutationList, UsherMetadata), Report>;
}

/// Convert graph to Usher MAT protobuf
pub fn usher_from_graph<N, E, D>(graph: &Graph<N, E, D>) -> Result<UsherTree, Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  (): UsherFromGraph<N, E, D>,
{
  let mut node_mutations = vec![];
  let mut condensed_nodes = vec![];
  let mut metadata = vec![];
  graph.iter_depth_first_preorder_forward(|node| {
    let edge = node.parents.first().map(|(_, edge)| edge.read_arc());
    let edge = edge.as_deref();
    let node = &node.payload;
    let (node, mutations, meta) =
      <() as UsherFromGraph<N, E, D>>::usher_node_from_graph_components(&UsherGraphContext { node, edge, graph })
        .unwrap();
    node_mutations.push(mutations);
    condensed_nodes.push(node);
    metadata.push(meta);
  });

  let newick = nwk_write_str(graph, &NwkWriteOptions::default())?;
  Ok(UsherTree {
    newick,
    node_mutations,
    condensed_nodes,
    metadata,
  })
}
