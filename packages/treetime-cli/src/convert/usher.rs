use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use eyre::Report;
use treetime::graph::graph::Graph;
use treetime::io::usher_mat::{
  UsherGraphContext, UsherMetadata, UsherMutationList, UsherRead, UsherTree, UsherTreeContext, UsherTreeNode,
  UsherWrite,
};
use treetime::make_internal_report;

pub struct UsherWriter;

impl UsherWrite<ConverterNode, ConverterEdge, ConverterData> for UsherWriter {
  fn new(graph: &Graph<ConverterNode, ConverterEdge, ConverterData>) -> Result<Self, Report> {
    Ok(Self {})
  }

  fn usher_node_from_graph_components(
    &mut self,
    context: &UsherGraphContext<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<(UsherTreeNode, UsherMutationList, UsherMetadata), Report> {
    let &UsherGraphContext { node, .. } = context;

    let node_name = node
      .name
      .as_ref()
      .ok_or_else(|| make_internal_report!("Encountered node with empty name"))?
      .to_owned();

    let usher_node = UsherTreeNode {
      node_name,
      condensed_leaves: vec![],
    };

    let usher_mutations = UsherMutationList { mutation: vec![] };

    let usher_meta = UsherMetadata {
      clade_annotations: vec![],
    };

    Ok((usher_node, usher_mutations, usher_meta))
  }
}

pub struct UsherReader;

impl UsherRead<ConverterNode, ConverterEdge, ConverterData> for UsherReader {
  fn new(tree: &UsherTree) -> Result<Self, Report> {
    Ok(Self {})
  }

  fn usher_data_to_graph_data(&mut self, tree: &UsherTree) -> Result<ConverterData, Report> {
    Ok(ConverterData::default())
  }

  fn usher_node_to_graph_components(
    &mut self,
    context: &UsherTreeContext,
  ) -> Result<(ConverterNode, ConverterEdge), Report> {
    let UsherTreeContext { node, .. } = context;

    let node = ConverterNode {
      name: node.name.as_ref().map(ToOwned::to_owned),
    };

    let edge = ConverterEdge { weight: None };

    Ok((node, edge))
  }
}
