#[cfg(test)]
mod tests {
  use crate::jobs::start_job;
  use app_commands::job::{JobEvent, JobId, JobRegistry, TerminalEvent};
  use helpers::{status, timetree_config_json};
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use std::sync::Arc;
  use std::thread;
  use tempfile::tempdir;
  use treetime_utils::assert_error;

  #[test]
  fn test_jobs_emit_started_first_and_end_with_one_terminal() {
    let registry = Arc::new(JobRegistry::default());
    let out = tempdir().unwrap();
    let job = start_job(&registry, "job-a", "timetree", &timetree_config_json(out.path())).unwrap();
    let events = Mutex::new(vec![]);
    let terminal = job.run(|event: JobEvent| events.lock().push(event));
    let events = events.into_inner();
    assert_eq!(
      (true, true, "ok", 0),
      (
        matches!(events.first(), Some(JobEvent::Started(started)) if started.job_id.as_str() == "job-a"),
        events.iter().any(|event| matches!(event, JobEvent::Progress(_))),
        status(&terminal),
        events
          .iter()
          .filter(|event| matches!(event, JobEvent::Terminal(_)))
          .count(),
      )
    );
    assert!(
      !registry.cancel(&JobId::parse("job-a").unwrap()),
      "a finished job leaves the registry"
    );
  }

  #[test]
  fn test_jobs_cancel_one_of_two_concurrent_jobs() {
    let registry = Arc::new(JobRegistry::default());
    let out_first = tempdir().unwrap();
    let out_second = tempdir().unwrap();
    let first = start_job(&registry, "first", "timetree", &timetree_config_json(out_first.path())).unwrap();
    let second = start_job(
      &registry,
      "second",
      "timetree",
      &timetree_config_json(out_second.path()),
    )
    .unwrap();

    let (first_terminal, second_terminal) = thread::scope(|scope| {
      let first_registry = Arc::clone(&registry);
      let first = scope.spawn(move || {
        first.run(move |event| {
          if matches!(event, JobEvent::Started(_)) {
            first_registry.cancel(&JobId::parse("first").unwrap());
          }
        })
      });
      let second = scope.spawn(move || second.run(|_event| {}));
      (first.join().unwrap(), second.join().unwrap())
    });
    assert_eq!(("cancelled", "ok"), (status(&first_terminal), status(&second_terminal)));
  }

  #[test]
  fn test_jobs_cancel_before_start_ends_cancelled() {
    let registry = Arc::new(JobRegistry::default());
    let out = tempdir().unwrap();
    let job = start_job(&registry, "early", "timetree", &timetree_config_json(out.path())).unwrap();
    assert!(registry.cancel(&JobId::parse("early").unwrap()));
    assert_eq!("cancelled", status(&job.run(|_event| {})));
  }

  #[test]
  fn test_jobs_reject_unknown_command() {
    let registry = Arc::new(JobRegistry::default());
    assert_error!(
      start_job(&registry, "x", "homoplasy", "{}"),
      "When reading the command name `homoplasy`: Matching variant not found"
    );
  }

  #[test]
  fn test_jobs_reject_duplicate_running_job_id() {
    let registry = Arc::new(JobRegistry::default());
    let _running = start_job(&registry, "same", "clock", "{}").unwrap();
    assert_error!(
      start_job(&registry, "same", "clock", "{}"),
      "a job with id `same` is already running"
    );
  }

  #[test]
  fn test_jobs_reject_malformed_config_json() {
    let registry = Arc::new(JobRegistry::default());
    assert_error!(
      start_job(&registry, "x", "clock", "{"),
      "When reading the command configuration: EOF while parsing an object at line 1 column 1"
    );
  }

  #[test]
  fn test_jobs_unknown_setting_ends_with_error_terminal() {
    let registry = Arc::new(JobRegistry::default());
    let job = start_job(&registry, "x", "clock", r#"{ "bogus": 1 }"#).unwrap();
    let terminal = job.run(|_event| {});
    let TerminalEvent::Error { message, .. } = terminal else {
      panic!("expected an error, got {terminal:?}");
    };
    assert_eq!("invalid configuration: unknown field `bogus`", message);
  }

  mod helpers {
    use app_commands::job::TerminalEvent;
    use serde_json::json;
    use std::path::Path;

    pub(super) fn status(terminal: &TerminalEvent) -> &'static str {
      match terminal {
        TerminalEvent::Ok { .. } => "ok",
        TerminalEvent::Error { .. } => "error",
        TerminalEvent::Cancelled { .. } => "cancelled",
      }
    }

    pub(super) fn timetree_config_json(out: &Path) -> String {
      let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      json!({
        "tree": zika.join("tree.nwk"),
        "metadata": zika.join("metadata.tsv"),
        "alignment": [zika.join("aln.fasta.xz")],
        "max_iter": 2,
        "seed": 7,
        "output_all": out,
      })
      .to_string()
    }
  }
}
