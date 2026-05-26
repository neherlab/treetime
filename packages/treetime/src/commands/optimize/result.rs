use crate::payload::ancestral::GraphAncestral;
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct OptimizeResult {
  #[serde(skip)]
  pub graph: GraphAncestral,
}
