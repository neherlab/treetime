use crate::graph::edge::GraphEdge;
use crate::graph::graph::Graph;
use crate::graph::node::{GraphNode, Named};
use crate::io::file::create_file_or_stdout;
use crate::io::file::open_file_or_stdin;
use crate::io::json::{
  json_read, json_read_file, json_read_str, json_write, json_write_file, json_write_str, JsonPretty,
};
use crate::io::nwk::{nwk_read_str, nwk_write_str, EdgeFromNwk, EdgeToNwk, NodeFromNwk, NodeToNwk, NwkWriteOptions};
use crate::make_error;
use bytes::Buf;
use eyre::{Report, WrapErr};
use smart_default::SmartDefault;
use std::io::{Read, Write};
use std::path::Path;

pub use usher_mat_utils::{UsherMetadata, UsherMutation, UsherMutationList, UsherTree, UsherTreeNode};

pub fn usher_mat_pb_read_file<C, N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let filepath = filepath.as_ref();
  usher_mat_pb_read::<C, _, _, _>(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading Usher MAT protobuf file '{filepath:#?}'"))
}

pub fn usher_mat_pb_read_bytes<C, N, E, D>(buf: impl Buf) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let tree = usher_mat_utils::usher_mat_pb_read_bytes(buf).wrap_err("When reading Usher MAT protobuf bytes")?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_pb_read<C, N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let tree = usher_mat_utils::usher_mat_pb_read(reader).wrap_err("When reading Usher MAT protobuf")?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_json_read_file<C, N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let filepath = filepath.as_ref();
  let tree = json_read_file(filepath).wrap_err_with(|| format!("When reading Usher MAT JSON file '{filepath:#?}'"))?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_json_read_str<C, N, E, D>(s: impl AsRef<str>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let tree = json_read_str(s).wrap_err("When reading Usher MAT JSON string")?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_json_read<C, N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode + NodeFromNwk + Named,
  E: GraphEdge + EdgeFromNwk,
  D: Sync + Send + Default,
  C: UsherRead<N, E, D>,
{
  let tree = json_read(reader).wrap_err("When reading Usher MAT JSON")?;
  usher_to_graph::<C, _, _, _>(&tree)
}

pub fn usher_mat_pb_write_file<C, N, E, D>(filepath: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  C: UsherWrite<N, E, D>,
{
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  usher_mat_pb_write::<C, _, _, _>(&mut f, graph)
    .wrap_err_with(|| format!("When reading Usher MAT protobuf file '{filepath:#?}'"))?;
  writeln!(f)?;
  Ok(())
}

pub fn usher_mat_pb_write_bytes<C, N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  C: UsherWrite<N, E, D>,
{
  let tree = usher_from_graph::<C, _, _, _>(graph)?;
  let mut buf = Vec::new();
  usher_mat_utils::usher_mat_pb_write_bytes(&mut buf, &tree).wrap_err("When writing Usher MAT protobuf bytes")
}

pub fn usher_mat_pb_write<C, N, E, D>(writer: &mut impl Write, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  C: UsherWrite<N, E, D>,
{
  let tree = usher_from_graph::<C, _, _, _>(graph)?;
  usher_mat_utils::usher_mat_pb_write(writer, &tree).wrap_err("When writing Usher MAT protobuf")
}

pub fn usher_mat_json_write_file<C, N, E, D>(
  filepath: impl AsRef<Path>,
  graph: &Graph<N, E, D>,

  options: &UsherMatJsonOptions,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  C: UsherWrite<N, E, D>,
{
  let filepath = filepath.as_ref();
  let tree = usher_from_graph::<C, _, _, _>(graph)?;
  json_write_file(filepath, &tree, JsonPretty(options.pretty))
    .wrap_err_with(|| format!("When writing Usher MAT JSON file: '{filepath:#?}'"))?;
  Ok(())
}

pub fn usher_mat_json_write_str<C, N, E, D>(
  graph: &Graph<N, E, D>,

  options: &UsherMatJsonOptions,
) -> Result<String, Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  C: UsherWrite<N, E, D>,
{
  let tree = usher_from_graph::<C, _, _, _>(graph)?;
  json_write_str(&tree, JsonPretty(options.pretty)).wrap_err("When writing Usher MAT JSON string")
}

pub fn usher_mat_json_write<C, N, E, D>(
  writer: &mut impl Write,
  graph: &Graph<N, E, D>,

  options: &UsherMatJsonOptions,
) -> Result<(), Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  C: UsherWrite<N, E, D>,
{
  let tree = usher_from_graph::<C, _, _, _>(graph)?;
  json_write(writer, &tree, JsonPretty(options.pretty)).wrap_err("When writing Usher MAT JSON")
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
  E: GraphEdge + EdgeFromNwk,
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
  E: GraphEdge + EdgeFromNwk,
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
    let context = UsherTreeContext {
      node: UsherNodeImpl {
        index: i,
        name: node.payload.name().map(|name| name.as_ref().to_owned()),
        branch_length: 0.0,
        clade_annotations: tree.metadata[i].clade_annotations.clone(),
        mutations: tree.node_mutations[i].mutation.clone(),
      },
      tree,
    };

    let (graph_node, _) = converter.usher_node_to_graph_components(&context).unwrap();
    *node.payload = graph_node;
    i += 1;
  });

  Ok(graph)
}

pub struct UsherGraphContext<'a, N, E, D>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  pub node: &'a N,
  pub edge: Option<&'a E>,
  pub graph: &'a Graph<N, E, D>,
}

pub trait UsherWrite<N, E, D>: Sized
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
{
  fn new(graph: &Graph<N, E, D>) -> Result<Self, Report>;

  fn usher_node_from_graph_components(
    &mut self,
    context: &UsherGraphContext<N, E, D>,
  ) -> Result<(UsherTreeNode, UsherMutationList, UsherMetadata), Report>;
}

/// Convert graph to Usher MAT protobuf
pub fn usher_from_graph<C, N, E, D>(graph: &Graph<N, E, D>) -> Result<UsherTree, Report>
where
  N: GraphNode + NodeToNwk,
  E: GraphEdge + EdgeToNwk,
  D: Sync + Send + Default,
  C: UsherWrite<N, E, D>,
{
  let mut node_mutations = vec![];
  let mut condensed_nodes = vec![];
  let mut metadata = vec![];

  let mut converter = C::new(graph)?;
  graph.iter_depth_first_preorder_forward(|node| {
    let edge = node.parents.first().map(|(_, edge)| edge.read_arc());
    let edge = edge.as_deref();
    let node = &node.payload;

    let (node, mutations, meta) = converter
      .usher_node_from_graph_components(&UsherGraphContext { node, edge, graph })
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
