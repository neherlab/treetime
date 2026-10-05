#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use helpers::{RecordingProgress, config_for};
  use itertools::Itertools;
  use pretty_assertions::assert_eq;
  use rstest::rstest;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;
  use treetime::progress::NoopProgress;
  use treetime_utils::pretty_assert_ulps_eq;

  #[rstest]
  #[case::timetree(AppCommand::Timetree)]
  #[case::clock(AppCommand::Clock)]
  #[case::ancestral_marginal(AppCommand::Ancestral)]
  #[case::mugration(AppCommand::Mugration)]
  #[case::optimize(AppCommand::Optimize)]
  #[case::prune(AppCommand::Prune)]
  #[trace]
  fn test_progress_fractions_increase_and_done_comes_last(#[case] command: AppCommand) {
    let outdir = tempdir().unwrap();
    let prepared = command.prepare_value(&config_for(command, outdir.path())).unwrap();
    let progress = RecordingProgress::new(outdir.path().to_path_buf());
    let outcome = prepared.args.run(&NoopCancel, &progress, &NoopProgress).unwrap();

    let events = progress.events();
    let fractions: Vec<f64> = events.iter().map(|event| event.fraction).collect();
    assert!(
      fractions
        .iter()
        .tuple_windows()
        .all(|(earlier, later)| earlier <= later),
      "stage fractions must not decrease: {events:?}"
    );

    let done: Vec<usize> = events
      .iter()
      .enumerate()
      .filter(|(_, event)| event.stage == "Done")
      .map(|(position, _)| position)
      .collect();
    assert_eq!(
      vec![events.len() - 1],
      done,
      "exactly one Done, as the last stage: {events:?}"
    );
    let last = events.last().unwrap();
    pretty_assert_ulps_eq!(1.0, last.fraction);

    let files_at_done = last.files.clone();
    assert!(!outcome.output_files.is_empty());
    assert_eq!(
      outcome
        .output_files
        .iter()
        .map(|file| file.path.clone())
        .collect::<Vec<_>>(),
      files_at_done,
      "every output is written before Done"
    );
  }

  #[test]
  fn test_progress_ancestral_parsimony_reports_done_once() {
    let outdir = tempdir().unwrap();
    let mut config = config_for(AppCommand::Ancestral, outdir.path());
    config["method_anc"] = "parsimony".into();
    let prepared = AppCommand::Ancestral.prepare_value(&config).unwrap();
    let progress = RecordingProgress::new(outdir.path().to_path_buf());
    prepared.args.run(&NoopCancel, &progress, &NoopProgress).unwrap();
    let stages: Vec<String> = progress.events().into_iter().map(|event| event.stage).collect();
    assert_eq!(
      vec![
        "Reading input",
        "Parsing tree",
        "Fitch parsimony",
        "Writing output",
        "Done"
      ],
      stages
    );
  }

  mod helpers {
    use crate::command::AppCommand;
    use parking_lot::Mutex;
    use serde_json::{Value, json};
    use std::fs;
    use std::path::{Path, PathBuf};
    use treetime::progress::StageSink;

    pub(super) fn config_for(command: AppCommand, outdir: &Path) -> Value {
      let data = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data");
      let zika = data.join("zika/20");
      let flu = data.join("flu/h3n2/20");
      match command {
        AppCommand::Timetree => json!({
          "tree": zika.join("tree.nwk"),
          "metadata": zika.join("metadata.tsv"),
          "alignment": [zika.join("aln.fasta.xz")],
          "max_iter": 2,
          "seed": 7,
          "output_all": outdir,
        }),
        AppCommand::Clock => json!({
          "tree": zika.join("tree.nwk"),
          "metadata": zika.join("metadata.tsv"),
          "output_all": outdir,
        }),
        AppCommand::Ancestral => json!({
          "tree": zika.join("tree.nwk"),
          "alignment": [zika.join("aln.fasta.xz")],
          "output_all": outdir,
        }),
        AppCommand::Mugration => json!({
          "tree": zika.join("tree.nwk"),
          "metadata": zika.join("metadata.tsv"),
          "attribute": "country",
          "output_all": outdir,
        }),
        AppCommand::Optimize => json!({
          "tree": flu.join("tree.nwk"),
          "alignment": [flu.join("aln.fasta.xz")],
          "output_all": outdir,
        }),
        AppCommand::Prune => json!({
          "tree": zika.join("tree.nwk"),
          "alignment": [zika.join("aln.fasta.xz")],
          "prune_short": 1e-6,
          "output_all": outdir,
        }),
      }
    }

    #[derive(Clone, Debug)]
    pub(super) struct RecordedStage {
      pub stage: String,
      pub fraction: f64,
      pub files: Vec<PathBuf>,
    }

    pub(super) struct RecordingProgress {
      outdir: PathBuf,
      events: Mutex<Vec<RecordedStage>>,
    }

    impl RecordingProgress {
      pub(super) fn new(outdir: PathBuf) -> Self {
        Self {
          outdir,
          events: Mutex::new(vec![]),
        }
      }

      pub(super) fn events(&self) -> Vec<RecordedStage> {
        self.events.lock().clone()
      }

      fn files(&self) -> Vec<PathBuf> {
        let mut files: Vec<PathBuf> = fs::read_dir(&self.outdir)
          .unwrap()
          .map(|entry| entry.unwrap().path())
          .filter(|path| path.is_file())
          .collect();
        files.sort();
        files
      }
    }

    impl StageSink for RecordingProgress {
      fn report(&self, stage: &str, fraction: f64, _message: &str) {
        let files = if stage == "Done" { self.files() } else { vec![] };
        self.events.lock().push(RecordedStage {
          stage: stage.to_owned(),
          fraction,
          files,
        });
      }
    }
  }
}
