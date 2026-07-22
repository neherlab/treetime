use crate::ancestral::pipeline::AncestralPartition;
use crate::commands::ancestral::aa_node_data::AaNodeData;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::payload::ancestral::GraphAncestral;
use serde::Serialize;

#[derive(Serialize)]
pub struct AncestralGraphData {
  pub partition: Option<AncestralPartition>,
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
  pub mask: Vec<bool>,
  pub aa_node_data: Option<AaNodeData>,
}

impl AncestralGraphData {
  pub fn new(
    partition: Option<AncestralPartition>,
    gtr: Option<GTR>,
    model_name: GtrModelName,
    mask: Vec<bool>,
    aa_node_data: Option<AaNodeData>,
  ) -> Self {
    Self {
      partition,
      gtr,
      model_name,
      mask,
      aa_node_data,
    }
  }
}

#[derive(Serialize)]
pub struct AncestralResult {
  #[serde(skip)]
  pub graph: GraphAncestral<AncestralGraphData>,
}

impl std::ops::Deref for AncestralResult {
  type Target = AncestralGraphData;

  fn deref(&self) -> &Self::Target {
    self.graph.data()
  }
}
