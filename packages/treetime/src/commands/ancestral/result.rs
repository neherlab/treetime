use crate::gtr::get_gtr::GtrModelName;
use crate::gtr::gtr::GTR;
use crate::payload::ancestral::GraphAncestral;

#[derive(Debug)]
pub struct AncestralResult {
  pub graph: GraphAncestral,
  pub gtr: Option<GTR>,
  pub model_name: GtrModelName,
}
