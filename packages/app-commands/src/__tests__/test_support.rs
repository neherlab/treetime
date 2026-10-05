#[cfg(test)]
pub(crate) mod tests {
  use crate::json_value::SparseConfig;
  use serde_json::Value;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack, Sub};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_primitives::AsciiChar;

  pub(crate) fn sparse(value: Value) -> SparseConfig {
    serde_json::from_value(value).expect("a test config is a mapping of settings")
  }

  pub(crate) fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .map(PathBuf::from)
      .expect("project has workspace root")
  }

  pub(crate) fn c(b: u8) -> AsciiChar {
    AsciiChar::from_byte_unchecked(b)
  }

  pub(crate) fn edge_mutation_map(
    graph: &Graph,
    edge_subs: &[(usize, Vec<Sub>)],
    edge_deletions: &[(usize, Vec<AlignedMutation>)],
  ) -> BTreeMap<GraphEdgeKey, Vec<Mutation>> {
    graph
      .get_edges()
      .enumerate()
      .map(|(idx, edge)| {
        let subs = edge_subs
          .iter()
          .filter(|(edge_idx, _)| *edge_idx == idx)
          .flat_map(|(_, subs)| subs.iter().cloned())
          .map(|sub| Mutation::substitution(MutationTrack::Nucleotide, sub));
        let deletions = edge_deletions
          .iter()
          .filter(|(edge_idx, _)| *edge_idx == idx)
          .flat_map(|(_, deletions)| deletions.iter())
          .map(|deletion| Mutation {
            track: MutationTrack::Nucleotide,
            event: MutationEvent::Deletion(deletion.clone()),
          });
        (edge.key(), subs.chain(deletions).collect())
      })
      .collect()
  }
}
