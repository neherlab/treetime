use crate::graph::edge::{Edge, GraphEdge};
use crate::graph::graph::Graph;
use crate::graph::node::GraphNode;
use crate::io::file::create_file_or_stdout;
use crate::io::file::open_file_or_stdin;
use crate::io::json::{
  JsonPretty, json_read, json_read_file, json_read_str, json_write, json_write_file, json_write_str,
};
use crate::make_internal_error;
use eyre::{Report, WrapErr};
use maplit::{btreemap, btreeset};
use parking_lot::RwLock;
pub use phyloxml::{
  PhyloXmlId, Phyloxml, PhyloxmlAccession, PhyloxmlAnnotation, PhyloxmlBinaryCharacterList, PhyloxmlBinaryCharacters,
  PhyloxmlBranchColor, PhyloxmlClade, PhyloxmlCladeRelation, PhyloxmlConfidence, PhyloxmlDate, PhyloxmlDistribution,
  PhyloxmlDomainArchitecture, PhyloxmlEvents, PhyloxmlMolSeq, PhyloxmlPhylogeny, PhyloxmlPoint, PhyloxmlProperty,
  PhyloxmlProteinDomain, PhyloxmlReference, PhyloxmlSequence, PhyloxmlSequenceRelation, PhyloxmlTaxonomy, PhyloxmlUri,
  Polygon,
};
use smart_default::SmartDefault;
use std::collections::VecDeque;
use std::io::{Cursor, Read, Write};
use std::path::Path;
use std::sync::Arc;

pub fn phyloxml_read_file<N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  (): PhyloxmlToGraph<N, E, D>,
{
  let filepath = filepath.as_ref();
  phyloxml_read(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading PhyloXML file '{}'", filepath.display()))
}

pub fn phyloxml_read_str<N, E, D>(s: impl AsRef<str>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  (): PhyloxmlToGraph<N, E, D>,
{
  phyloxml_read(Cursor::new(s.as_ref())).wrap_err("When reading PhyloXML string")
}

pub fn phyloxml_read<N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  (): PhyloxmlToGraph<N, E, D>,
{
  let pxml = phyloxml::phyloxml_read(reader).wrap_err("When reading PhyloXML")?;
  phyloxml_to_graph(&pxml)
}

pub fn phyloxml_json_read_file<N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  (): PhyloxmlToGraph<N, E, D>,
{
  let filepath = filepath.as_ref();
  let tree =
    json_read_file(filepath).wrap_err_with(|| format!("When reading PhyloXML JSON file '{}'", filepath.display()))?;
  phyloxml_to_graph(&tree)
}

pub fn phyloxml_json_read_str<N, E, D>(s: impl AsRef<str>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  (): PhyloxmlToGraph<N, E, D>,
{
  let pxml = json_read_str(s).wrap_err("When reading PhyloXML JSON string")?;
  phyloxml_to_graph(&pxml)
}

pub fn phyloxml_json_read<N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  (): PhyloxmlToGraph<N, E, D>,
{
  let pxml = json_read(reader).wrap_err("When reading PhyloXML JSON")?;
  phyloxml_to_graph(&pxml)
}

pub fn phyloxml_write_file<N, E, D>(filepath: impl AsRef<Path>, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataFromGraphData + Sync + Send,
  (): PhyloxmlFromGraph<N, E, D>,
{
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  phyloxml_write(&mut f, graph).wrap_err_with(|| format!("When reading PhyloXML file '{}'", filepath.display()))?;
  writeln!(f)?;
  Ok(())
}

pub fn phyloxml_write_str<N, E, D>(graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataFromGraphData + Sync + Send,
  (): PhyloxmlFromGraph<N, E, D>,
{
  let mut buf = Vec::new();
  phyloxml_write(&mut buf, graph).wrap_err("When writing PhyloXML string")
}

pub fn phyloxml_write<N, E, D>(writer: &mut impl Write, graph: &Graph<N, E, D>) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataFromGraphData + Sync + Send,
  (): PhyloxmlFromGraph<N, E, D>,
{
  let pxml = phyloxml_from_graph(graph)?;
  phyloxml::phyloxml_write(writer, &pxml).wrap_err("When writing PhyloXML")
}

pub fn phyloxml_json_write_file<N, E, D>(
  filepath: impl AsRef<Path>,
  graph: &Graph<N, E, D>,
  options: &PhyloxmlJsonOptions,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataFromGraphData + Sync + Send,
  (): PhyloxmlFromGraph<N, E, D>,
{
  let filepath = filepath.as_ref();
  let pxml = phyloxml_from_graph(graph)?;
  json_write_file(filepath, &pxml, JsonPretty(options.pretty))
    .wrap_err_with(|| format!("When writing PhyloXML JSON file: '{}'", filepath.display()))?;
  Ok(())
}

pub fn phyloxml_json_write_str<N, E, D>(graph: &Graph<N, E, D>, options: &PhyloxmlJsonOptions) -> Result<String, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataFromGraphData + Sync + Send,
  (): PhyloxmlFromGraph<N, E, D>,
{
  let pxml = phyloxml_from_graph(graph)?;
  json_write_str(&pxml, JsonPretty(options.pretty)).wrap_err("When writing PhyloXML JSON string")
}

pub fn phyloxml_json_write<N, E, D>(
  writer: &mut impl Write,
  graph: &Graph<N, E, D>,
  options: &PhyloxmlJsonOptions,
) -> Result<(), Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataFromGraphData + Sync + Send,
  (): PhyloxmlFromGraph<N, E, D>,
{
  let pxml = phyloxml_from_graph(graph)?;
  json_write(writer, &pxml, JsonPretty(options.pretty)).wrap_err("When writing PhyloXML JSON")
}

pub trait PhyloxmlDataToGraphData: Sized {
  fn phyloxml_data_to_graph_data(pxml: &Phyloxml) -> Result<Self, Report>;
}

pub trait PhyloxmlDataFromGraphData {
  fn phyloxml_data_from_graph_data(&self) -> Result<Phyloxml, Report>;
}

#[derive(SmartDefault)]
pub struct PhyloxmlJsonOptions {
  #[default = true]
  pretty: bool,
}

pub struct PhyloxmlContext<'a> {
  pub clade: &'a PhyloxmlClade,
  pub pxml: &'a Phyloxml,
}

/// Describes conversion from Phyloxml tree node data when reading from PhyloXML
pub trait PhyloxmlToGraph<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
{
  fn phyloxml_node_to_graph_components(context: &PhyloxmlContext) -> Result<(N, E), Report>;
}

pub struct PhyloxmlNodeImpl {
  pub name: Option<String>,
  pub branch_length: f64,
}

/// Convert PhyloXML to graph
pub fn phyloxml_to_graph<N, E, D>(pxml: &Phyloxml) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  (): PhyloxmlToGraph<N, E, D>,
{
  let phylogeny = {
    let n_phylogeny = pxml.phylogeny.len();
    if n_phylogeny != 1 {
      make_internal_error!("Only PhyloXML with exactly one phylogeny are currently supported, but found {n_phylogeny}")
    } else {
      Ok(&pxml.phylogeny[0])
    }
  }?;

  let mut graph = Graph::<N, E, D>::with_data(D::phyloxml_data_to_graph_data(pxml)?);
  let mut queue: VecDeque<_> = phylogeny.clade.iter().map(|clade| (None, clade)).collect();
  while let Some((parent_key, clade)) = queue.pop_front() {
    let (graph_node, graph_edge) =
      <() as PhyloxmlToGraph<N, E, D>>::phyloxml_node_to_graph_components(&PhyloxmlContext { clade, pxml })?;
    let node_key = graph.add_node(graph_node);
    if let Some(parent_key) = parent_key {
      graph.add_edge(parent_key, node_key, graph_edge)?;
    }
    for child in &clade.clade {
      queue.push_back((Some(node_key), child));
    }
  }
  graph.build()?;

  Ok(graph)
}

pub struct PhyloxmlGraphContext<'a, N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataFromGraphData + Sync + Send,
  (): PhyloxmlFromGraph<N, E, D>,
{
  pub node: &'a N,
  pub edge: Option<&'a E>,
  pub graph: &'a Graph<N, E, D>,
}

/// Describes conversion to Phyloxml tree node data when writing to PhyloXML
pub trait PhyloxmlFromGraph<N, E, D>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataFromGraphData + Sync + Send,
  (): PhyloxmlFromGraph<N, E, D>,
{
  fn phyloxml_node_from_graph_components(context: &PhyloxmlGraphContext<N, E, D>) -> Result<PhyloxmlClade, Report>;
}

/// Convert graph to PhyloXML
pub fn phyloxml_from_graph<N, E, D>(graph: &Graph<N, E, D>) -> Result<Phyloxml, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataFromGraphData + Sync + Send,
  (): PhyloxmlFromGraph<N, E, D>,
{
  let root = graph.get_exactly_one_root()?;

  // Pre-order iteration to construct the nodes from graph nodes and edges
  let mut node_map = {
    let mut node_map = btreemap! {};
    let mut queue = VecDeque::from([(Arc::clone(&root), None)]);
    while let Some((current_node, current_edge)) = queue.pop_front() {
      let current_node = current_node.read_arc();

      let node = &*current_node.payload().read_arc();
      let edge = current_edge
        .as_ref()
        .map(|edge: &Arc<RwLock<Edge<E>>>| edge.read_arc().payload().read_arc());
      let edge = edge.as_deref();
      let current_tree_node =
        <() as PhyloxmlFromGraph<N, E, D>>::phyloxml_node_from_graph_components(&PhyloxmlGraphContext {
          node,
          edge,
          graph,
        })?;
      for (child, edge) in graph.children_of(&current_node) {
        queue.push_back((child, Some(edge)));
      }
      node_map.insert(current_node.key(), current_tree_node);
    }
    node_map
  };

  // Post-order traversal to populate .children array
  let mut visited = btreeset! {};
  let mut stack = vec![Arc::clone(&root)];
  while let Some(node) = stack.pop() {
    if visited.contains(&node.read_arc().key()) {
      let mut children_to_add = vec![];
      for (child, _) in graph.children_of(&node.read_arc()) {
        let child_key = child.read_arc().key();
        if let Some(child_node) = node_map.remove(&child_key) {
          children_to_add.push(child_node);
        }
      }
      if let Some(node) = node_map.get_mut(&node.read_arc().key()) {
        node.clade.extend(children_to_add);
      }
    } else {
      visited.insert(node.read_arc().key());
      stack.push(Arc::clone(&node));
      for (child, _) in graph.children_of(&node.read_arc()) {
        stack.push(child);
      }
    }
  }

  let mut data = graph.data().read_arc().phyloxml_data_from_graph_data()?;
  let clade = node_map.remove(&root.read_arc().key()).unwrap();
  data.phylogeny[0].clade = Some(clade);
  Ok(data)
}
