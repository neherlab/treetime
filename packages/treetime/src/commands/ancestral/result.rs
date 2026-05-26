use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::payload::ancestral::GraphAncestral;
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct AncestralResult {
  #[serde(skip)]
  pub graph: GraphAncestral,
  #[serde(skip)]
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
}
