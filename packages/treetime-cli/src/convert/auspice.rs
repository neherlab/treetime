use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use crate::convert::mutation::{PartitionedMutations, format_mutation_list, parse_mutation_list};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use serde_json::Value;
use std::collections::BTreeMap;
use treetime::graph::graph::Graph;
use treetime::graph::node::GraphNodeKey;
use treetime::io::auspice::{AuspiceGraphContext, AuspiceRead, AuspiceTreeContext, AuspiceWrite};
use treetime::io::auspice_types::{
  AuspiceTree, AuspiceTreeBranchAttrs, AuspiceTreeData, AuspiceTreeNode, AuspiceTreeNodeAttrs,
};
use treetime::make_internal_report;

fn parse_auspice_mutations(auspice_mutations: &BTreeMap<String, Vec<String>>) -> Result<PartitionedMutations, Report> {
  auspice_mutations
    .iter()
    .map(|(partition, muts)| Ok((partition.clone(), parse_mutation_list(muts)?)))
    .try_collect()
}

fn format_auspice_mutations(mutations: &PartitionedMutations) -> BTreeMap<String, Vec<String>> {
  mutations
    .iter()
    .map(|(partition, muts)| (partition.clone(), format_mutation_list(muts)))
    .collect()
}

pub struct AuspiceWriter {
  pub divs: BTreeMap<GraphNodeKey, f64>,
}

impl AuspiceWrite<ConverterNode, ConverterEdge, ConverterData> for AuspiceWriter {
  fn new(_: &Graph<ConverterNode, ConverterEdge, ConverterData>) -> Result<Self, Report> {
    Ok(Self { divs: btreemap! {} })
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
    let &AuspiceGraphContext {
      node_key,
      node,
      parent_key,
      parent,
      edge,
      graph,
    } = context;

    let name = node
      .name
      .as_ref()
      .ok_or_else(|| make_internal_report!("Encountered node with empty name"))?
      .to_owned();

    let parent_div = parent_key
      .and_then(|parent_key| self.divs.get(&parent_key).copied())
      .or(Some(0.0)); // FIXME: we probably should not assume 0 if parent div is missing or if it's a root node

    let branch_length = edge.and_then(|edge| edge.branch_length);

    let div = if let (Some(parent_div), Some(branch_length)) = (parent_div, branch_length) {
      Some(parent_div + branch_length)
    } else {
      Some(0.0) // FIXME: we probably should not assume 0 if div or branch_length is missing
    };

    if let Some(div) = div {
      self.divs.insert(node_key, div);
    }

    let mutations = edge.map_or(BTreeMap::new(), |edge| format_auspice_mutations(&edge.mutations));

    Ok(AuspiceTreeNode {
      name,
      branch_attrs: AuspiceTreeBranchAttrs {
        mutations,
        labels: None,
        other: Value::default(),
      },
      node_attrs: AuspiceTreeNodeAttrs {
        div,
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

pub struct AuspiceReader {
  has_mutations: bool,
}

impl AuspiceRead<ConverterNode, ConverterEdge, ConverterData> for AuspiceReader {
  fn new(tree: &AuspiceTree) -> Result<Self, Report> {
    let has_mutations = tree
      .iter_depth_first_preorder()
      .any(|(_, node)| !node.branch_attrs.mutations.is_empty());
    Ok(Self { has_mutations })
  }

  fn auspice_data_to_graph_data(&mut self, tree: &AuspiceTree) -> Result<ConverterData, Report> {
    Ok(ConverterData {
      rooted: true,
      has_mutations: self.has_mutations,
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
    let branch_length = context.branch_length();
    let mutations = parse_auspice_mutations(&node.branch_attrs.mutations)?;
    Ok((
      ConverterNode {
        name: Some(node.name.clone()),
      },
      ConverterEdge {
        branch_length,
        mutations,
      },
    ))
  }
}
