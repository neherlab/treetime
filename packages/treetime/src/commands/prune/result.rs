use crate::payload::ancestral::GraphAncestral;
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct PruneResult {
  #[serde(skip)]
  pub graph: GraphAncestral,
}
