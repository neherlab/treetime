use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use eyre::Report;
use serde_json::Value;
use treetime::graph::graph::Graph;
use treetime::io::auspice::{
  AuspiceGraphContext, AuspiceRead, AuspiceTree, AuspiceTreeBranchAttrs, AuspiceTreeContext, AuspiceTreeData,
  AuspiceTreeNode, AuspiceTreeNodeAttrs, AuspiceWrite,
};

pub struct AuspiceWriter {}

impl AuspiceWrite<ConverterNode, ConverterEdge, ConverterData> for AuspiceWriter {
  fn new(graph: &Graph<ConverterNode, ConverterEdge, ConverterData>) -> Result<Self, Report> {
    Ok(Self {})
  }

  fn auspice_data_from_graph_data(
    &self,
    graph: &Graph<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<AuspiceTreeData, Report> {
    let data = graph.data().read_arc();
    Ok(AuspiceTreeData {
      version: data.version.clone(),
      meta: data.meta.clone(),
      root_sequence: data.root_sequence.clone(),
      other: data.other.clone(),
    })
  }

  fn auspice_node_from_graph_components(
    &mut self,
    context: &AuspiceGraphContext<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<AuspiceTreeNode, Report> {
    let AuspiceGraphContext { node, edge, .. } = context;
    Ok(AuspiceTreeNode {
      name: node.name.clone().unwrap_or_default(),
      branch_attrs: AuspiceTreeBranchAttrs::default(),
      node_attrs: AuspiceTreeNodeAttrs {
        div: edge.and_then(|edge| edge.weight),
        clade_membership: None,
        region: None,
        country: None,
        division: None,
        other: Value::default(),
      },
      children: vec![],
      other: Value::default(),
    })
  }
}

pub struct AuspiceReader {}

impl AuspiceRead<ConverterNode, ConverterEdge, ConverterData> for AuspiceReader {
  fn new(tree: &AuspiceTree) -> Result<Self, Report> {
    Ok(Self {})
  }

  fn auspice_data_to_graph_data(&mut self, tree: &AuspiceTree) -> Result<ConverterData, Report> {
    Ok(ConverterData {
      rooted: true,
      version: tree.data.version.clone(),
      meta: tree.data.meta.clone(),
      root_sequence: tree.data.root_sequence.clone(),
      other: tree.data.other.clone(),
    })
  }

  fn auspice_node_to_graph_components(
    &mut self,
    context: &AuspiceTreeContext,
  ) -> Result<(ConverterNode, ConverterEdge), Report> {
    let AuspiceTreeContext { node, .. } = context;
    Ok((
      ConverterNode {
        name: Some(node.name.clone()),
      },
      ConverterEdge {
        weight: node.node_attrs.div,
      },
    ))
  }
}
