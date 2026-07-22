use crate::nwk::{EdgeFromNwk, NodeFromNwk, nwk_read_str};
use bytes::{Buf, BytesMut};
use eyre::{Report, WrapErr};
use smart_default::SmartDefault;
use std::io::{Read, Write};
use std::path::Path;
use treetime_graph::edge::{GraphEdge, HasBranchLength};
use treetime_graph::graph::Graph;
use treetime_graph::node::{GraphNode, Named};
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::io::json::{
  JsonPretty, json_read, json_read_file, json_read_str, json_write, json_write_file, json_write_str,
};
use treetime_utils::make_error;

pub use util_usher_mat::{UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode};

pub fn usher_mat_pb_read_file<C, N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk + HasBranchLength,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let filepath = filepath.as_ref();
  usher_mat_pb_read::<C, _, _, _>(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading Usher MAT protobuf file '{}'", filepath.display()))
}

pub fn usher_mat_pb_read_bytes<C, N, E, D>(buf: impl Buf) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk + HasBranchLength,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let tree = util_usher_mat::usher_mat_pb_read_bytes(buf).wrap_err("When reading Usher MAT protobuf bytes")?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_pb_read<C, N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk + HasBranchLength,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let tree = util_usher_mat::usher_mat_pb_read(reader).wrap_err("When reading Usher MAT protobuf")?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_json_read_file<C, N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk + HasBranchLength,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let filepath = filepath.as_ref();
  let tree =
    json_read_file(filepath).wrap_err_with(|| format!("When reading Usher MAT JSON file '{}'", filepath.display()))?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_json_read_str<C, N, E, D>(s: impl AsRef<str>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk + HasBranchLength,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let tree = json_read_str(s).wrap_err("When reading Usher MAT JSON string")?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_json_read<C, N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk + HasBranchLength,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let tree = json_read(reader).wrap_err("When reading Usher MAT JSON")?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_pb_write_file(filepath: impl AsRef<Path>, tree: &UsherTree) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  usher_mat_pb_write(&mut f, tree)
    .wrap_err_with(|| format!("When writing Usher MAT protobuf file '{}'", filepath.display()))?;
  writeln!(f)?;
  Ok(())
}

pub fn usher_mat_pb_write_bytes(tree: &UsherTree) -> Result<Vec<u8>, Report> {
  let mut buf = BytesMut::new();
  util_usher_mat::usher_mat_pb_write_bytes(&mut buf, tree).wrap_err("When writing Usher MAT protobuf bytes")?;
  Ok(buf.to_vec())
}

pub fn usher_mat_pb_write(writer: &mut impl Write, tree: &UsherTree) -> Result<(), Report> {
  util_usher_mat::usher_mat_pb_write(writer, tree).wrap_err("When writing Usher MAT protobuf")
}

pub fn usher_mat_json_write_file(
  filepath: impl AsRef<Path>,
  tree: &UsherTree,
  options: &UsherMatJsonOptions,
) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  json_write_file(filepath, tree, JsonPretty(options.pretty))
    .wrap_err_with(|| format!("When writing Usher MAT JSON file: '{}'", filepath.display()))?;
  Ok(())
}

pub fn usher_mat_json_write_str(tree: &UsherTree, options: &UsherMatJsonOptions) -> Result<String, Report> {
  json_write_str(tree, JsonPretty(options.pretty)).wrap_err("When writing Usher MAT JSON string")
}

pub fn usher_mat_json_write(
  writer: &mut impl Write,
  tree: &UsherTree,
  options: &UsherMatJsonOptions,
) -> Result<(), Report> {
  json_write(writer, tree, JsonPretty(options.pretty)).wrap_err("When writing Usher MAT JSON")
}

#[derive(SmartDefault)]
pub struct UsherMatJsonOptions {
  #[default = true]
  pretty: bool,
}

pub struct UsherNodeImpl {
  pub index: usize,
  pub name: Option<String>,
  pub branch_length: f64,
  pub clade_annotations: Vec<String>,
  pub mutations: Vec<UsherMutation>,
}

pub struct UsherTreeContext<'a> {
  pub node: UsherNodeImpl,
  pub tree: &'a UsherTree,
}

pub trait UsherRead<N, E, D>: Sized
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk + HasBranchLength,
  D: Sync + Send + Default,
{
  fn new(tree: &UsherTree) -> Result<Self, Report>;

  fn usher_data_to_graph_data(&mut self, tree: &UsherTree) -> Result<D, Report>;

  fn usher_node_to_graph_components(&mut self, context: &UsherTreeContext) -> Result<(N, E), Report>;
}

/// Convert Usher MAT protobuf to graph
pub fn usher_to_graph<C, N, E, D>(tree: &UsherTree) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk + HasBranchLength,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let mut graph: Graph<N, E, D> = nwk_read_str(&tree.newick)?;

  let n_nodes = &graph.num_nodes();
  let n_muts = &tree.node_mutations.len();
  let n_meta = &tree.metadata.len();

  if n_nodes != n_muts || n_nodes != n_meta {
    return make_error!(
      "Inconsistent number of nodes, mutations and metadata entries: nodes: {n_nodes}, muts: {n_muts}, meta: {n_meta}"
    );
  }

  let mut converter = C::new(tree)?;

  graph.set_data(converter.usher_data_to_graph_data(tree)?);

  let mut i = 0;
  graph.iter_depth_first_preorder_forward(|mut node| {
    // Branch lengths come from the embedded Newick parsed above; recover them so the
    // converter can place them on the reconstructed edge alongside the mutations.
    let branch_length = node
      .parents
      .first()
      .and_then(|(_, edge)| edge.read_arc().branch_length())
      .unwrap_or(0.0);
    let context = UsherTreeContext {
      node: UsherNodeImpl {
        index: i,
        name: node.payload.name().map(|name| name.as_ref().to_owned()),
        branch_length,
        clade_annotations: tree.metadata[i].clade_annotations.clone(),
        mutations: tree.node_mutations[i].mutation.clone(),
      },
      tree,
    };

    let (graph_node, graph_edge) = converter.usher_node_to_graph_components(&context)?;
    *node.payload = graph_node;
    if let Some((_, edge)) = node.parents.first() {
      *edge.write_arc() = graph_edge;
    }
    i += 1;
    Ok(())
  })?;

  Ok(graph)
}
