#[cfg(test)]
pub(crate) mod tests {
  use eyre::Report;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use treetime::seq::indel::InDel;
  use treetime::seq::mutation::{Mutation, MutationTrack, Sub};
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
    edge_indels: &[(usize, Vec<InDel>)],
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
        let indels = edge_indels
          .iter()
          .filter(|(edge_idx, _)| *edge_idx == idx)
          .flat_map(|(_, indels)| indels.iter())
          .map(|indel| Mutation::indel(MutationTrack::Nucleotide, indel));
        Ok((edge.key(), subs.chain(indels).collect::<Result<Vec<_>, Report>>()?))
      })
      .collect()
  }
}
