use crate::payload::ancestral::GraphAncestral;

#[derive(Debug)]
pub struct PruneResult {
  pub graph: GraphAncestral,
}
