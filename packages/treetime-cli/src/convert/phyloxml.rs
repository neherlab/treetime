use crate::convert::convert::{ConverterData, ConverterEdge, ConverterNode};
use eyre::Report;
use serde_json::Value;
use std::collections::BTreeMap;
use treetime::io::auspice::AuspiceTreeMeta;
use treetime::io::phyloxml::{
  Phyloxml, PhyloxmlClade, PhyloxmlContext, PhyloxmlDataFromGraphData, PhyloxmlDataToGraphData, PhyloxmlFromGraph,
  PhyloxmlGraphContext, PhyloxmlPhylogeny, PhyloxmlToGraph,
};

impl PhyloxmlDataFromGraphData for ConverterData {
  fn phyloxml_data_from_graph_data(&self) -> Result<Phyloxml, Report> {
    Ok(Phyloxml {
      phylogeny: vec![PhyloxmlPhylogeny {
        rooted: self.rooted,
        rerootable: None,
        branch_length_unit: None,
        phylogeny_type: None,
        name: None,
        id: None,
        description: None,
        date: None,
        confidence: vec![],
        clade: None,
        clade_relation: vec![],
        sequence_relation: vec![],
        property: vec![],
        other: BTreeMap::default(),
      }],
      other: BTreeMap::default(),
    })
  }
}

impl PhyloxmlDataToGraphData for ConverterData {
  fn phyloxml_data_to_graph_data(pxml: &Phyloxml) -> Result<Self, Report> {
    Ok(Self {
      rooted: pxml.phylogeny[0].rooted,
      version: None,
      meta: AuspiceTreeMeta::default(),
      root_sequence: None,
      other: Value::default(),
    })
  }
}

impl PhyloxmlFromGraph<ConverterNode, ConverterEdge, ConverterData> for () {
  fn phyloxml_node_from_graph_components(
    PhyloxmlGraphContext { edge, .. }: &PhyloxmlGraphContext<ConverterNode, ConverterEdge, ConverterData>,
  ) -> Result<PhyloxmlClade, Report> {
    Ok(PhyloxmlClade {
      name: None,
      branch_length_elem: edge.and_then(|edge| edge.weight),
      branch_length_attr: edge.and_then(|edge| edge.weight),
      confidence: vec![],
      width: None,
      color: None,
      node_id: None,
      taxonomy: vec![],
      sequence: vec![],
      events: None,
      binary_characters: None,
      distribution: vec![],
      date: None,
      reference: vec![],
      property: vec![],
      clade: vec![],
      other: BTreeMap::default(),
    })
  }
}

impl PhyloxmlToGraph<ConverterNode, ConverterEdge, ConverterData> for () {
  fn phyloxml_node_to_graph_components(
    PhyloxmlContext { clade, pxml }: &PhyloxmlContext,
  ) -> Result<(ConverterNode, ConverterEdge), Report> {
    Ok((
      ConverterNode {
        name: clade
          .name
          .clone()
          .or_else(|| clade.taxonomy.first().as_ref().and_then(|t| t.scientific_name.clone())),
      },
      ConverterEdge {
        weight: clade.branch_length_attr.or(clade.branch_length_elem),
      },
    ))
  }
}
