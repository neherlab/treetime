#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::job::{CancelToken, JobEvent, JobId, JobProgress, TerminalEvent, run_job};
  use eyre::Report;
  use helpers::{status, timetree_config};
  use parking_lot::Mutex;
  use pretty_assertions::{assert_eq, assert_ne};
  use rstest::rstest;
  use serde_json::json;
  use std::thread;
  use tempfile::tempdir;
  use treetime::cancel::{Cancel, NoopCancel};
  use treetime::progress::{LogSink, NoopProgress, StageSink};
  use treetime_utils::assert_error;

  #[test]
  fn test_job_success_ends_ok_with_output_files() {
    let outdir = tempdir().unwrap();
    let job_id = JobId::parse("job-1").unwrap();
    let terminal = run_job(&job_id, &NoopCancel, || {
      AppCommand::Timetree
        .prepare_value(&timetree_config(outdir.path()))?
        .args
        .run(&NoopCancel, &NoopProgress, &NoopProgress)
    });
    let TerminalEvent::Ok { job_id: id, result } = terminal else {
      panic!("expected ok, got {terminal:?}");
    };
    assert_eq!(
      (job_id, AppCommand::Timetree, true),
      (
        id,
        result.command,
        result
          .output_files
          .iter()
          .any(|file| file.path == outdir.path().join("timetree.auspice.json"))
      )
    );
  }

  #[test]
  fn test_job_rejected_config_ends_error_with_cli_message() {
    let job_id = JobId::parse("job-2").unwrap();
    let terminal = run_job(&job_id, &NoopCancel, || {
      AppCommand::Clock
        .prepare_value(&json!({ "tree": "t.nwk", "no_such_setting": 1 }))?
        .args
        .run(&NoopCancel, &NoopProgress, &NoopProgress)
    });
    assert_eq!(
      json!({
        "status": "error",
        "job_id": "job-2",
        "message": "invalid configuration: unknown field `no_such_setting`",
        "causes": [],
      }),
      serde_json::to_value(terminal).unwrap()
    );
  }

  #[test]
  fn test_job_command_failure_ends_error_with_cause_chain() {
    let outdir = tempdir().unwrap();
    let mut config = timetree_config(outdir.path());
    config["metadata"] = json!(outdir.path().join("missing.tsv"));
    let terminal = run_job(&JobId::parse("job-3").unwrap(), &NoopCancel, || {
      AppCommand::Timetree
        .prepare_value(&config)?
        .args
        .run(&NoopCancel, &NoopProgress, &NoopProgress)
    });
    let TerminalEvent::Error { message, causes, .. } = terminal else {
      panic!("expected error, got {terminal:?}");
    };
    assert!(!message.is_empty());
    assert!(
      causes.iter().any(|cause| cause.contains("missing.tsv")),
      "the cause chain names the missing file: {message:?} {causes:?}"
    );
  }

  #[test]
  fn test_job_cancelled_before_start_ends_cancelled() {
    let outdir = tempdir().unwrap();
    let token = CancelToken::default();
    token.cancel();
    let terminal = run_job(&JobId::parse("job-4").unwrap(), &token, || {
      token.check()?;
      AppCommand::Timetree
        .prepare_value(&timetree_config(outdir.path()))?
        .args
        .run(&token, &NoopProgress, &NoopProgress)
    });
    assert_eq!(
      json!({ "status": "cancelled", "job_id": "job-4" }),
      serde_json::to_value(terminal).unwrap()
    );
  }

  #[test]
  fn test_job_panic_ends_error() {
    let terminal = run_job(&JobId::parse("job-5").unwrap(), &NoopCancel, || panic!("boom"));
    assert_eq!(
      json!({
        "status": "error",
        "job_id": "job-5",
        "message": "internal error: the computation panicked: boom",
        "causes": [],
      }),
      serde_json::to_value(terminal).unwrap()
    );
  }

  #[test]
  fn test_job_error_ends_error_with_message() {
    let terminal = run_job(&JobId::parse("job-6").unwrap(), &NoopCancel, || {
      Err(Report::msg("input path is outside the examples folder"))
    });
    assert_eq!(
      json!({
        "status": "error",
        "job_id": "job-6",
        "message": "input path is outside the examples folder",
        "causes": [],
      }),
      serde_json::to_value(terminal).unwrap()
    );
  }

  #[test]
  fn test_job_cancelling_one_of_two_concurrent_jobs_leaves_the_other_running() {
    let first = CancelToken::default();
    let second = CancelToken::default();
    let dir_first = tempdir().unwrap();
    let dir_second = tempdir().unwrap();

    let (first_terminal, second_terminal) = thread::scope(|scope| {
      let first_job = scope.spawn(|| {
        let cancel_on_first_event = JobProgress::new(|_event| first.cancel());
        run_job(&JobId::parse("first").unwrap(), &first, || {
          AppCommand::Timetree
            .prepare_value(&timetree_config(dir_first.path()))?
            .args
            .run(&first, &cancel_on_first_event, &cancel_on_first_event)
        })
      });
      let second_job = scope.spawn(|| {
        run_job(&JobId::parse("second").unwrap(), &second, || {
          AppCommand::Timetree
            .prepare_value(&timetree_config(dir_second.path()))?
            .args
            .run(&second, &NoopProgress, &NoopProgress)
        })
      });
      (first_job.join().unwrap(), second_job.join().unwrap())
    });

    assert_eq!(
      ("cancelled", "ok", true, false),
      (
        status(&first_terminal),
        status(&second_terminal),
        first.is_cancelled(),
        second.is_cancelled()
      )
    );
  }

  #[test]
  fn test_job_id_rejects_unsafe_characters() {
    assert_error!(
      JobId::parse("../x"),
      "invalid job id `../x`: expected 1 to 128 ASCII letters, digits, `-` or `_`"
    );
  }

  #[rustfmt::skip]
  #[rstest]
  #[case::traversal(r#""../../etc""#, "invalid job id `../../etc`: expected 1 to 128 ASCII letters, digits, `-` or `_`")]
  #[case::slash(    r#""a/b""#,       "invalid job id `a/b`: expected 1 to 128 ASCII letters, digits, `-` or `_`")]
  #[case::empty(    r#""""#,          "invalid job id ``: expected 1 to 128 ASCII letters, digits, `-` or `_`")]
  #[trace]
  fn test_job_id_deserialization_rejects_invalid_ids(#[case] json: &str, #[case] expected: &str) {
    assert_error!(serde_json::from_str::<JobId>(json).map_err(Report::new), expected);
  }

  #[test]
  fn test_job_id_deserialization_accepts_a_valid_id() {
    let id: JobId = serde_json::from_str(r#""run_1-a""#).unwrap();
    assert_eq!("run_1-a", id.as_str());
  }

  #[test]
  fn test_job_random_ids_differ() {
    assert_ne!(JobId::random(), JobId::random());
  }

  #[test]
  fn test_job_progress_forwards_progress_and_log_events() {
    let events = Mutex::new(vec![]);
    let progress = JobProgress::new(|event: JobEvent| events.lock().push(serde_json::to_value(event).unwrap()));
    progress.report("Reading input", 0.0, "");
    treetime::progress_warn!(progress, "tip {} has no date", "A");
    assert_eq!(
      vec![
        json!({ "type": "progress", "data": { "stage": "Reading input", "fraction": 0.0, "message": "" } }),
        json!({ "type": "log", "data": { "level": "warn", "message": "tip A has no date" } }),
      ],
      *events.lock()
    );
  }

  mod helpers {
    use crate::job::TerminalEvent;
    use serde_json::{Value, json};
    use std::path::Path;

    pub(super) fn status(terminal: &TerminalEvent) -> &'static str {
      match terminal {
        TerminalEvent::Ok { .. } => "ok",
        TerminalEvent::Error { .. } => "error",
        TerminalEvent::Cancelled { .. } => "cancelled",
        TerminalEvent::Interrupted { .. } => "interrupted",
      }
    }

    pub(super) fn timetree_config(outdir: &Path) -> Value {
      let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      json!({
        "tree": zika.join("tree.nwk"),
        "metadata": zika.join("metadata.tsv"),
        "alignment": [zika.join("aln.fasta.xz")],
        "max_iter": 2,
        "seed": 7,
        "output_all": outdir,
      })
    }
  }
}
