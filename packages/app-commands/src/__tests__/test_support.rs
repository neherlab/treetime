#[cfg(test)]
pub(crate) mod tests {
  use deser::Deserialize;
  use eyre::Report;
  use itertools::izip;
  use ndarray::{Array1, Array2, array};
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use std::collections::BTreeMap;
  use std::path::{Path, PathBuf};
  use std::sync::LazyLock;
  use treetime::alphabet::alphabet::{Alphabet, AlphabetName};
  use treetime::gtr::get_gtr::{GtrModel, GtrModelName, HKY85Params};
  use treetime::gtr::gtr::GTR;
  use treetime::progress::{LogLevel, LogSink, RunWarning};
  use treetime::seq::mutation::{AlignedMutation, Mutation, MutationEvent, MutationTrack, Sub};
  use treetime_graph::edge::GraphEdgeKey;
  use treetime_graph::graph::Graph;
  use treetime_io::fasta::{FastaWriter, fasta_read_file};
  use treetime_primitives::{AsciiChar, Seq};
  use treetime_utils::adapters::{Array2Rows, ArrayVec};
  use treetime_utils::pretty_assert_ulps_eq;

  pub(crate) const HKY85_PARAMS: [&str; 2] = ["kappa=0.2", "pis=0.1,0.2,0.3,0.4"];

  pub(crate) static NUC_ALPHABET: LazyLock<Alphabet> = LazyLock::new(Alphabet::default);

  pub(crate) fn project_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
      .parent()
      .and_then(|p| p.parent())
      .map(PathBuf::from)
      .expect("project has workspace root")
  }

  pub(crate) fn write_column_blocks(alignment: &Path, cut: usize, dir: &Path) -> Result<[PathBuf; 2], Report> {
    let records = fasta_read_file(alignment, &*NUC_ALPHABET)?;
    let paths = [dir.join("columns-a.fasta"), dir.join("columns-b.fasta")];
    for (path, first_block) in izip!(&paths, [true, false]) {
      let mut writer = FastaWriter::create(path)?;
      for record in &records {
        let (start, end) = if first_block { (0, cut) } else { (cut, record.seq.len()) };
        writer.write(
          &record.seq_name,
          record.desc.as_deref(),
          &Seq::from(&record.seq[start..end]),
        )?;
      }
      writer.finish()?;
    }
    Ok(paths)
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

  pub(crate) fn hky85_model() -> GtrModel {
    GtrModel::HKY85(HKY85Params {
      kappa: 0.2,
      pi: Some(array![0.1, 0.2, 0.3, 0.4]),
      ..HKY85Params::default()
    })
  }

  pub(crate) fn assert_written_gtr(expected: &GtrModel, actual: &WrittenGtr) -> Result<(), Report> {
    let gtr: GTR = expected.build(AlphabetName::Nuc)?;
    assert_eq!(expected.name(), actual.model_name);
    pretty_assert_ulps_eq!(gtr.pi().clone(), actual.pi.clone(), max_ulps = 4);
    pretty_assert_ulps_eq!(gtr.W().clone(), actual.w.clone(), max_ulps = 4);
    Ok(())
  }

  #[derive(Debug, Deserialize)]
  pub(crate) struct WrittenModels {
    pub(crate) models: Vec<WrittenGtr>,
  }

  #[derive(Debug, Deserialize)]
  pub(crate) struct WrittenGtr {
    pub(crate) model_name: GtrModelName,
    #[deser(deserialize_as = ArrayVec)]
    pub(crate) pi: Array1<f64>,
    #[deser(rename = "W", deserialize_as = Array2Rows)]
    pub(crate) w: Array2<f64>,
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
