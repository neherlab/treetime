#[cfg(test)]
pub(crate) mod tests {
  use eyre::Report;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack, Sub};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_primitives::AsciiChar;

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
  ) -> Result<BTreeMap<GraphEdgeKey, Vec<Mutation>>, Report> {
    graph
      .get_edges()
      .enumerate()
      .map(|(idx, edge)| {
        let subs = edge_subs
          .iter()
          .filter(|(edge_idx, _)| *edge_idx == idx)
          .flat_map(|(_, subs)| subs.iter().cloned())
          .map(|sub| Ok(Mutation::substitution(MutationTrack::Nucleotide, sub)));
        let deletions = edge_deletions
          .iter()
          .filter(|(edge_idx, _)| *edge_idx == idx)
          .flat_map(|(_, deletions)| deletions.iter())
          .map(|deletion| {
            Ok(Mutation {
              track: MutationTrack::Nucleotide,
              event: MutationEvent::Deletion(deletion.clone()),
            })
          });
        Ok((edge.key(), subs.chain(deletions).collect::<Result<Vec<_>, Report>>()?))
      })
      .collect()
  }
}
