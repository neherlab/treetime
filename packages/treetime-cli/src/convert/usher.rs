use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use crate::convert::mutation::{Mutation, MutationList, PartitionedMutations};
use eyre::Report;
use itertools::Itertools;
use maplit::btreemap;
use treetime::make_internal_report;
use treetime_graph::graph::Graph;
use treetime_io::usher_mat::{
  UsherGraphContext, UsherMetadata, UsherMutation, UsherMutationList, UsherRead, UsherTree, UsherTreeContext,
  UsherTreeNode, UsherWrite,
};

const NUC_PARTITION: &str = "nuc";
const NUC_CHARS: [u8; 4] = [b'A', b'C', b'G', b'T'];

fn nuc_int_to_char(nuc: i32) -> Option<u8> {
  NUC_CHARS.get(nuc as usize).copied()
}

fn nuc_char_to_int(nuc: u8) -> Option<i32> {
  NUC_CHARS.iter().position(|&c| c == nuc).map(|i| i as i32)
}

fn parse_usher_mutations(usher_mutations: &[UsherMutation]) -> Result<PartitionedMutations, Report> {
  if usher_mutations.is_empty() {
    return Ok(btreemap! {});
  }

  let mutations: MutationList = usher_mutations
    .iter()
    .filter_map(|m| {
      let reference = nuc_int_to_char(m.par_nuc)?;
      let alternative = m.mut_nuc.first().and_then(|&n| nuc_int_to_char(n))?;
      let position = m.position as usize;
      Some(Mutation::new(reference, position, alternative))
    })
    .collect_vec();

  if mutations.is_empty() {
    Ok(btreemap! {})
  } else {
    Ok(btreemap! { NUC_PARTITION.to_owned() => mutations })
  }
}

fn format_usher_mutations(mutations: &PartitionedMutations) -> UsherMutationList {
  let nuc_mutations = mutations.get(NUC_PARTITION);

  let mutation = nuc_mutations.map_or(vec![], |muts| {
    muts
      .iter()
      .filter_map(|m| {
        let par_nuc = nuc_char_to_int(m.reference)?;
        let mut_nuc = nuc_char_to_int(m.alternative)?;
        Some(UsherMutation {
          position: m.position as i32,
          ref_nuc: par_nuc,
          par_nuc,
          mut_nuc: vec![mut_nuc],
          chromosome: String::new(),
        })
      })
      .collect_vec()
  });

  UsherMutationList { mutation }
}

pub struct UsherWriter;

impl UsherWrite<ConverterNode, ConverterEdge, ConverterData> for UsherWriter {
  fn new(_graph: &Graph<ConverterNode, ConverterEdge, ConverterData>) -> Result<Self, Report> {
    Ok(Self {})
  }

  fn usher_node_from_graph_components(
    &mut self,
    context: &UsherGraphContext<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<(UsherTreeNode, UsherMutationList, UsherMetadata), Report> {
    let &UsherGraphContext { node, edge, .. } = context;

    let node_name = node
      .name
      .as_ref()
      .ok_or_else(|| make_internal_report!("Encountered node with empty name"))?
      .to_owned();

    let usher_node = UsherTreeNode {
      node_name,
      condensed_leaves: vec![],
    };

    let usher_mutations = edge.map_or(UsherMutationList { mutation: vec![] }, |e| {
      format_usher_mutations(&e.mutations)
    });

    let usher_meta = UsherMetadata {
      clade_annotations: vec![],
    };

    Ok((usher_node, usher_mutations, usher_meta))
  }
}

pub struct UsherReader {
  has_mutations: bool,
}

impl UsherRead<ConverterNode, ConverterEdge, ConverterData> for UsherReader {
  fn new(tree: &UsherTree) -> Result<Self, Report> {
    let has_mutations = tree.node_mutations.iter().any(|ml| !ml.mutation.is_empty());
    Ok(Self { has_mutations })
  }

  fn usher_data_to_graph_data(&mut self, _tree: &UsherTree) -> Result<ConverterData, Report> {
    Ok(ConverterData {
      has_mutations: self.has_mutations,
      ..Default::default()
    })
  }

  fn usher_node_to_graph_components(
    &mut self,
    context: &UsherTreeContext,
  ) -> Result<(ConverterNode, ConverterEdge), Report> {
    let UsherTreeContext { node, .. } = context;

    let mutations = parse_usher_mutations(&node.mutations)?;

    let node = ConverterNode {
      name: node.name.as_ref().map(ToOwned::to_owned),
    };

    let edge = ConverterEdge {
      branch_length: None,
      mutations,
    };

    Ok((node, edge))
  }
}
