use crate::ancestral::pipeline::AncestralPartition;
use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::payload::ancestral::GraphAncestral;
use serde::Serialize;

#[allow(clippy::manual_non_exhaustive, clippy::partial_pub_fields)] // The private unit field preserves Graph JSON's `data: null` shape.
#[derive(Serialize)]
#[serde(transparent)]
pub struct AncestralGraphData {
  marker: (),
  #[serde(skip)]
  pub partition: Option<AncestralPartition>,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  #[serde(skip)]
  pub model_name: GtrModelName,
  #[serde(skip)]
  pub mask: Vec<bool>,
}

impl AncestralGraphData {
  pub fn new(
    partition: Option<AncestralPartition>,
    gtr: Option<GTR>,
    model_name: GtrModelName,
    mask: Vec<bool>,
  ) -> Self {
    Self {
      marker: (),
      partition,
      gtr,
      model_name,
      mask,
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
