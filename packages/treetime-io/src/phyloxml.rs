use eyre::{Report, WrapErr};
use smart_default::SmartDefault;
use std::collections::VecDeque;
use std::io::{Cursor, Read, Write};
use std::path::Path;
use treetime_graph::edge::GraphEdge;
use treetime_graph::graph::Graph;
use treetime_graph::node::GraphNode;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::file::open_file_or_stdin;
use treetime_utils::io::json::{
  JsonPretty, json_read, json_read_file, json_read_str, json_write, json_write_file, json_write_str,
};
use treetime_utils::make_internal_error;
pub use util_phyloxml::{
  Phyloxml, PhyloxmlAccession, PhyloxmlAnnotation, PhyloxmlBinaryCharacterList, PhyloxmlBinaryCharacters,
  PhyloxmlBranchColor, PhyloxmlClade, PhyloxmlCladeRelation, PhyloxmlConfidence, PhyloxmlDate, PhyloxmlDistribution,
  PhyloxmlDomainArchitecture, PhyloxmlEvents, PhyloxmlId, PhyloxmlMolSeq, PhyloxmlPhylogeny, PhyloxmlPoint,
  PhyloxmlProperty, PhyloxmlProteinDomain, PhyloxmlReference, PhyloxmlSequence, PhyloxmlSequenceRelation,
  PhyloxmlTaxonomy, PhyloxmlUri, Polygon,
};

pub fn phyloxml_read_file<C, N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  C: PhyloxmlToGraph<N, E, D>,
{
  let filepath = filepath.as_ref();
  phyloxml_read::<C, _, _, _>(open_file_or_stdin(&Some(filepath))?)
    .wrap_err_with(|| format!("When reading PhyloXML file '{}'", filepath.display()))
}

pub fn phyloxml_read_str<C, N, E, D>(s: impl AsRef<str>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  C: PhyloxmlToGraph<N, E, D>,
{
  phyloxml_read::<C, _, _, _>(Cursor::new(s.as_ref())).wrap_err("When reading PhyloXML string")
}

pub fn phyloxml_read<C, N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  C: PhyloxmlToGraph<N, E, D>,
{
  let pxml = util_phyloxml::phyloxml_read(reader).wrap_err("When reading PhyloXML")?;
  phyloxml_to_graph::<C, _, _, _>(&pxml)
}

pub fn phyloxml_json_read_file<C, N, E, D>(filepath: impl AsRef<Path>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  C: PhyloxmlToGraph<N, E, D>,
{
  let filepath = filepath.as_ref();
  let tree =
    json_read_file(filepath).wrap_err_with(|| format!("When reading PhyloXML JSON file '{}'", filepath.display()))?;
  phyloxml_to_graph::<C, _, _, _>(&tree)
}

pub fn phyloxml_json_read_str<C, N, E, D>(s: impl AsRef<str>) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  C: PhyloxmlToGraph<N, E, D>,
{
  let pxml = json_read_str(s).wrap_err("When reading PhyloXML JSON string")?;
  phyloxml_to_graph::<C, _, _, _>(&pxml)
}

pub fn phyloxml_json_read<C, N, E, D>(reader: impl Read) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  C: PhyloxmlToGraph<N, E, D>,
{
  let pxml = json_read(reader).wrap_err("When reading PhyloXML JSON")?;
  phyloxml_to_graph::<C, _, _, _>(&pxml)
}

pub fn phyloxml_write_file(filepath: impl AsRef<Path>, phyloxml: &Phyloxml) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  phyloxml_write(&mut f, phyloxml).wrap_err_with(|| format!("When writing PhyloXML file '{}'", filepath.display()))?;
  writeln!(f)?;
  Ok(())
}

pub fn phyloxml_write_str(phyloxml: &Phyloxml) -> Result<String, Report> {
  let mut buf = Vec::new();
  phyloxml_write(&mut buf, phyloxml).wrap_err("When writing PhyloXML string")?;
  String::from_utf8(buf).wrap_err("PhyloXML output is not valid UTF-8")
}

pub fn phyloxml_write(writer: &mut impl Write, phyloxml: &Phyloxml) -> Result<(), Report> {
  util_phyloxml::phyloxml_write(writer, phyloxml).wrap_err("When writing PhyloXML")
}

pub fn phyloxml_json_write_file(
  filepath: impl AsRef<Path>,
  phyloxml: &Phyloxml,
  options: &PhyloxmlJsonOptions,
) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  json_write_file(filepath, phyloxml, JsonPretty(options.pretty))
    .wrap_err_with(|| format!("When writing PhyloXML JSON file: '{}'", filepath.display()))?;
  Ok(())
}

pub fn phyloxml_json_write_str(phyloxml: &Phyloxml, options: &PhyloxmlJsonOptions) -> Result<String, Report> {
  json_write_str(phyloxml, JsonPretty(options.pretty)).wrap_err("When writing PhyloXML JSON string")
}

pub fn phyloxml_json_write(
  writer: &mut impl Write,
  phyloxml: &Phyloxml,
  options: &PhyloxmlJsonOptions,
) -> Result<(), Report> {
  json_write(writer, phyloxml, JsonPretty(options.pretty)).wrap_err("When writing PhyloXML JSON")
}

pub trait PhyloxmlDataToGraphData: Sized {
  fn phyloxml_data_to_graph_data(pxml: &Phyloxml) -> Result<Self, Report>;
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
pub fn phyloxml_to_graph<C, N, E, D>(pxml: &Phyloxml) -> Result<Graph<N, E, D>, Report>
where
  N: GraphNode,
  E: GraphEdge,
  D: PhyloxmlDataToGraphData + Sync + Send,
  C: PhyloxmlToGraph<N, E, D>,
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
    let (graph_node, graph_edge) = C::phyloxml_node_to_graph_components(&PhyloxmlContext { clade, pxml })?;
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
