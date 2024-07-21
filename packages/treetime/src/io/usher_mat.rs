use crate::graph::assign_node_names::assign_node_names;
use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::io::file::create_file_or_stdout;
use crate::io::file::open_file_or_stdin;
use bytes::Buf;
use eyre::{Report, WrapErr};
use std::io::{Read, Write};
use std::path::Path;
use usher_mat_utils::parsimony::{CondensedNode, Data};

pub fn usher_mat_pb_read_file<N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: FromUsherMatData + Sync + Send,
  (N, E): FromUsherMatNode,
{
  let filepath = filepath.as_ref();
  usher_mat_pb_read(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading Usher MAT protobuf file '{filepath:#?}'"))
}

pub fn usher_mat_pb_read_bytes<N, E, D>(buf: impl Buf) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: FromUsherMatData + Sync + Send,
  (N, E): FromUsherMatNode,
{
  let data = usher_mat_utils::usher_mat_pb_read_bytes(buf).wrap_err("When reading Usher MAT protobuf bytes")?;
  convert_usher_mat_to_graph(&data)
}

pub fn usher_mat_pb_read<N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + Named,
  E: GraphEdge,
  D: FromUsherMatData + Sync + Send,
  (N, E): FromUsherMatNode,
{
  let data = usher_mat_utils::usher_mat_pb_read(reader).wrap_err("When reading Usher MAT protobuf")?;
  convert_usher_mat_to_graph(&data)
}

pub fn usher_mat_pb_write_file<N, E, D>(filepath: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: ToUsherMatData + Sync + Send,
  for<'a> (&'a N, Option<&'a E>): ToUsherMatNode,
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
  N: GraphNode,
  E: GraphEdge,
  D: ToUsherMatData + Sync + Send,
  for<'a> (&'a N, Option<&'a E>): ToUsherMatNode,
{
  let data = convert_graph_to_usher_mat(graph)?;
  let mut buf = Vec::new();
  usher_mat_utils::usher_mat_pb_write_bytes(&mut buf, &data).wrap_err("When writing Usher MAT protobuf bytes")
}

pub fn usher_mat_pb_write<N, E, D>(writer: &mut impl Write, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: ToUsherMatData + Sync + Send,
  for<'a> (&'a N, Option<&'a E>): ToUsherMatNode,
{
  let data = convert_graph_to_usher_mat(graph)?;
  usher_mat_utils::usher_mat_pb_write(writer, &data).wrap_err("When writing Usher MAT protobuf")
}

/// Defines conversion to tree node when writing to Usher MAT protobuf
pub trait ToUsherMatNode {
  fn to_usher_node(&self) -> CondensedNode;
}

/// Defines conversion from tree node when reading from Usher MAT protobuf
pub trait FromUsherMatNode: Sized {
  fn from_usher_node(mat_node: &CondensedNode) -> Self;
}

/// Defines conversion to tree global data when writing to Usher MAT protobuf
pub trait ToUsherMatData {
  fn to_usher_data(&self) -> Data;
}

/// Defines conversion from tree global data when reading from Usher MAT protobuf
pub trait FromUsherMatData: Sized {
  fn from_usher_data(data: &Data) -> Self;
}

fn convert_usher_mat_to_graph<N, E, D>(data: &Data) -> Result<Graph<N, E, D>, Report>
where
  (N, E): FromUsherMatNode,
  D: FromUsherMatData + Send + Sync,
  E: GraphEdge,
  N: GraphNode + Named,
{
  let mut graph = Graph::<N, E, D>::with_data(D::from_usher_data(data));

  graph.build()?;
  assign_node_names(&graph);
  Ok(graph)
}

fn convert_graph_to_usher_mat<N, E, D>(graph: &Graph<N, E, D>) -> Result<Data, Report>
where
  D: Send + Sync + ToUsherMatData,
  E: GraphEdge,
  N: GraphNode,
  for<'a> (&'a N, Option<&'a E>): ToUsherMatNode,
{
  let newick = String::new();
  // let newick = nwk_write_str(graph, &NwkWriteOptions::default())?;

  Ok(Data {
    newick,
    node_mutations: vec![],
    condensed_nodes: vec![],
    metadata: vec![],
  })
}
