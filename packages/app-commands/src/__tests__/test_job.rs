#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::job::{CancelToken, JobEvent, JobId, JobProgress, JobRegistry, TerminalEvent, run_job};
  use eyre::Report;
  use helpers::{accept, status, timetree_config};
  use parking_lot::Mutex;
  use pretty_assertions::{assert_eq, assert_ne};
  use serde_json::json;
  use std::sync::Arc;
  use std::thread;
  use tempfile::tempdir;
  use treetime::cancel::{Cancel, NoopCancel};
  use treetime::progress::{NoopProgress, ProgressSink};
  use treetime_utils::assert_error;

  #[test]
  fn test_job_success_ends_ok_with_output_files() {
    let outdir = tempdir().unwrap();
    let job_id = JobId::parse("job-1").unwrap();
    let terminal = run_job(
      &job_id,
      AppCommand::Timetree,
      &timetree_config(outdir.path()),
      &accept,
      &NoopCancel,
      &NoopProgress,
    );
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
          .contains(&outdir.path().join("timetree.auspice.json"))
      )
    );
  }

  #[test]
  fn test_job_rejected_config_ends_error_with_cli_message() {
    let job_id = JobId::parse("job-2").unwrap();
    let terminal = run_job(
      &job_id,
      AppCommand::Clock,
      &json!({ "tree": "t.nwk", "no_such_setting": 1 }),
      &accept,
      &NoopCancel,
      &NoopProgress,
    );
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
    let terminal = run_job(
      &JobId::parse("job-3").unwrap(),
      AppCommand::Timetree,
      &config,
      &accept,
      &NoopCancel,
      &NoopProgress,
    );
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
    let terminal = run_job(
      &JobId::parse("job-4").unwrap(),
      AppCommand::Timetree,
      &timetree_config(outdir.path()),
      &accept,
      &token,
      &NoopProgress,
    );
    assert_eq!(
      json!({ "status": "cancelled", "job_id": "job-4" }),
      serde_json::to_value(terminal).unwrap()
    );
  }

  #[test]
  fn test_job_panic_ends_error() {
    let terminal = run_job(
      &JobId::parse("job-5").unwrap(),
      AppCommand::Timetree,
      &json!({}),
      &|_| panic!("boom"),
      &NoopCancel,
      &NoopProgress,
    );
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
  fn test_job_prepare_hook_error_ends_error() {
    let terminal = run_job(
      &JobId::parse("job-6").unwrap(),
      AppCommand::Timetree,
      &json!({}),
      &|_| Err(Report::msg("input path is outside the data directory")),
      &NoopCancel,
      &NoopProgress,
    );
    assert_eq!(
      json!({
        "status": "error",
        "job_id": "job-6",
        "message": "input path is outside the data directory",
        "causes": [],
      }),
      serde_json::to_value(terminal).unwrap()
    );
  }

  #[test]
  fn test_job_cancelling_one_of_two_concurrent_jobs_leaves_the_other_running() {
    let registry = Arc::new(JobRegistry::default());
    let first = registry.register(JobId::parse("first").unwrap()).unwrap();
    let second = registry.register(JobId::parse("second").unwrap()).unwrap();
    let dir_first = tempdir().unwrap();
    let dir_second = tempdir().unwrap();

    let (first_terminal, second_terminal) = thread::scope(|scope| {
      let registry_for_first = Arc::clone(&registry);
      let first_job = scope.spawn(|| {
        let cancel_on_first_event = JobProgress::new(move |_event| {
          registry_for_first.cancel(&JobId::parse("first").unwrap());
        });
        run_job(
          first.job_id(),
          AppCommand::Timetree,
          &timetree_config(dir_first.path()),
          &accept,
          first.token(),
          &cancel_on_first_event,
        )
      });
      let second_job = scope.spawn(|| {
        run_job(
          second.job_id(),
          AppCommand::Timetree,
          &timetree_config(dir_second.path()),
          &accept,
          second.token(),
          &NoopProgress,
        )
      });
      (first_job.join().unwrap(), second_job.join().unwrap())
    });

    assert_eq!(
      ("cancelled", "ok", true, false),
      (
        status(&first_terminal),
        status(&second_terminal),
        first.token().is_cancelled(),
        second.token().is_cancelled()
      )
    );
  }

  #[test]
  fn test_job_registry_rejects_duplicate_id_and_forgets_finished_jobs() {
    let registry = Arc::new(JobRegistry::default());
    let job_id = JobId::parse("dup").unwrap();
    let handle = registry.register(job_id.clone()).unwrap();
    assert_error!(
      registry.register(job_id.clone()),
      "a job with id `dup` is already running"
    );
    drop(handle);
    assert!(!registry.cancel(&job_id), "a finished job leaves the registry");
    registry.register(job_id).unwrap();
  }

  #[test]
  fn test_job_registry_cancel_is_idempotent_and_scoped_to_one_job() {
    let registry = Arc::new(JobRegistry::default());
    let first = registry.register(JobId::parse("a").unwrap()).unwrap();
    let second = registry.register(JobId::parse("b").unwrap()).unwrap();
    let requests = (registry.cancel(first.job_id()), registry.cancel(first.job_id()));
    let third = registry.register(JobId::parse("c").unwrap()).unwrap();
    assert_eq!(
      ((true, true), true, false, false),
      (
        requests,
        first.token().is_cancelled(),
        second.token().is_cancelled(),
        third.token().is_cancelled()
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
    use eyre::Report;
    use serde_json::{Value, json};
    use std::path::Path;

    pub(super) fn status(terminal: &TerminalEvent) -> &'static str {
      match terminal {
        TerminalEvent::Ok { .. } => "ok",
        TerminalEvent::Error { .. } => "error",
        TerminalEvent::Cancelled { .. } => "cancelled",
      }
    }

    pub(super) fn accept(_config: &mut Value) -> Result<(), Report> {
      Ok(())
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
