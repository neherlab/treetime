#[cfg(test)]
pub(crate) mod tests {
  use eyre::Report;
  use maplit::btreemap;
  use std::collections::BTreeMap;
  use std::iter;
  use std::path::PathBuf;
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::ancestral::pipeline::SparseReconstruction;
  use treetime::gtr::get_gtr::{JC69Params, jc69};
  use treetime::partition::marginal::shared::update::MarginalEdges;
  use treetime::partition::marginal::sparse::partition::PartitionMarginalSparse;
  use treetime::partition::storage::sparse::{SparseEdgeObs, SparseNodeObs, SparseNodeState};
  use treetime::seq::mutation::{Mutation, MutationTrack, Sub};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_primitives::{AsciiChar, Seq};

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
    partition: &SparseReconstruction,
  ) -> Result<BTreeMap<GraphEdgeKey, Vec<Mutation>>, Report> {
    graph
      .get_edges()
      .map(|edge| {
        let key = edge.key();
        Ok((key, partition.edge_mutations(key, &MutationTrack::Nucleotide)?))
      })
      .collect()
  }

  pub(crate) fn make_test_partition(
    graph: &Graph,
    length: usize,
    edge_subs: &[(usize, Vec<Sub>)],
  ) -> Result<SparseReconstruction, Report> {
    let alphabet = Alphabet::default();
    let mut ref_seq: Seq = iter::repeat_with(|| c(b'A')).take(length).collect();
    for (_, subs) in edge_subs {
      for s in subs {
        if s.pos() < length {
          ref_seq[s.pos()] = s.reff();
        }
      }
    }

    let mut obs_nodes = btreemap! {};
    let mut node_states = btreemap! {};
    for node in graph.get_nodes() {
      let key = node.key();
      obs_nodes.insert(key, SparseNodeObs::new(&ref_seq, &alphabet));
      node_states.insert(key, SparseNodeState::leaf(&ref_seq));
    }

    let mut obs_edges = btreemap! {};
    let mut estimates = btreemap! {};
    let edges = graph.get_edges().collect::<Vec<_>>();
    for (idx, subs) in edge_subs {
      if let Some(edge) = edges.get(*idx) {
        let edge_key = edge.key();
        obs_edges.insert(edge_key, SparseEdgeObs::with_fitch_subs(subs.clone()));
        estimates.insert(edge_key, subs.clone());
      }
    }

    let partition = PartitionMarginalSparse {
      index: 0,
      alphabet,
      length,
      root_sequence: ref_seq,
      obs_nodes,
      obs_edges,
    };

    Ok(SparseReconstruction {
      partition,
      gtr: jc69(JC69Params::default())?,
      node_states,
      edges: MarginalEdges {
        estimates,
        ..MarginalEdges::default()
      },
    })
  }
}
