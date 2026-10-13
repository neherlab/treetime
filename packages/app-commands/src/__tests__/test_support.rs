#[cfg(test)]
pub(crate) mod tests {
  use parking_lot::Mutex;
  use std::collections::BTreeMap;
  use std::path::PathBuf;
  use std::sync::LazyLock;
  use treetime::alphabet::alphabet::Alphabet;
  use treetime::progress::{LogLevel, LogSink, RunWarning};
  use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack, Sub};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_primitives::AsciiChar;

  pub(crate) static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

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

  #[derive(Default)]
  pub(crate) struct RecordingLog {
    events: Mutex<Vec<(LogLevel, String)>>,
    warnings: Mutex<Vec<RunWarning>>,
  }

  impl RecordingLog {
    pub(crate) fn events(&self) -> Vec<(LogLevel, String)> {
      self.events.lock().clone()
    }

    pub(crate) fn messages(&self, level: LogLevel) -> Vec<String> {
      self
        .events()
        .into_iter()
        .filter(|(event_level, _)| *event_level == level)
        .map(|(_, message)| message)
        .collect()
    }

    pub(crate) fn into_warnings(self) -> Vec<RunWarning> {
      self.warnings.into_inner()
    }
  }

  impl LogSink for RecordingLog {
    fn log(&self, level: LogLevel, message: &str) {
      self.events.lock().push((level, message.to_owned()));
    }

    fn log_enabled(&self, _level: LogLevel) -> bool {
      true
    }

    fn warning(&self, warning: &RunWarning) {
      self.warnings.lock().push(warning.clone());
      self.log(LogLevel::Warn, &warning.message);
    }
  }
}
