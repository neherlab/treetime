#[cfg(test)]
mod tests {
  use crate::command::AppCommand;
  use crate::job::{IterationEvent, JobEvent, JobProgress};
  use helpers::{timetree_config, tracelog_rows};
  use parking_lot::Mutex;
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use tempfile::tempdir;
  use treetime::cancel::NoopCancel;

  #[test]
  fn test_iteration_events_match_the_tracelog_rows() {
    let outdir = tempdir().unwrap();
    let prepared = AppCommand::Timetree
      .prepare_value(&timetree_config(outdir.path()))
      .unwrap();
    let events = Mutex::new(vec![]);
    let progress = JobProgress::new(|event| {
      if let JobEvent::Iteration(iteration) = event {
        events.lock().push(iteration);
      }
    });
    prepared.args.run(&NoopCancel, &progress).unwrap();
    let events = events.into_inner();

    let from_events: Vec<Vec<String>> = events
      .iter()
      .map(|event| {
        let value = serde_json::to_value(event).unwrap();
        [
          "n_diff",
          "n_resolved",
          "max_time_change",
          "rms_time_change",
          "log_lh_seq",
          "log_lh_pos",
          "log_lh_coal",
          "log_lh_total",
        ]
        .iter()
        .map(|key| match &value[key] {
          serde_json::Value::Null => String::new(),
          other => other.as_f64().unwrap().to_string(),
        })
        .collect()
      })
      .collect();
    assert_eq!(tracelog_rows(&outdir.path().join("timetree.tracelog.csv")), from_events);
    assert_eq!(
      (1..=events.len()).collect::<Vec<_>>(),
      events.iter().map(|event| event.iteration).collect::<Vec<_>>()
    );
  }

  #[test]
  fn test_iteration_events_carry_the_estimated_clock_rate_and_r_squared() {
    let outdir = tempdir().unwrap();
    let prepared = AppCommand::Timetree
      .prepare_value(&timetree_config(outdir.path()))
      .unwrap();
    let events: Mutex<Vec<IterationEvent>> = Mutex::new(vec![]);
    let progress = JobProgress::new(|event| {
      if let JobEvent::Iteration(iteration) = event {
        events.lock().push(iteration);
      }
    });
    prepared.args.run(&NoopCancel, &progress).unwrap();
    let events = events.into_inner();
    assert!(!events.is_empty());
    assert!(
      events.iter().all(|event| event.clock_rate.0 > 0.0
        && event
          .r_squared
          .is_some_and(|r_squared| (0.0..=1.0).contains(&r_squared.0))),
      "every iteration of an estimated clock has a positive rate and an R² in [0, 1]: {events:?}"
    );
  }

  #[test]
  fn test_iteration_events_of_a_fixed_rate_have_that_rate_and_no_r_squared() {
    let outdir = tempdir().unwrap();
    let mut config = timetree_config(outdir.path());
    config["clock_rate"] = json!(0.0008);
    let prepared = AppCommand::Timetree.prepare_value(&config).unwrap();
    let events: Mutex<Vec<IterationEvent>> = Mutex::new(vec![]);
    let progress = JobProgress::new(|event| {
      if let JobEvent::Iteration(iteration) = event {
        events.lock().push(iteration);
      }
    });
    prepared.args.run(&NoopCancel, &progress).unwrap();
    let summary: Vec<(f64, bool)> = events
      .into_inner()
      .iter()
      .map(|event| (event.clock_rate.0, event.r_squared.is_none()))
      .collect();
    assert_eq!(vec![(0.0008, true), (0.0008, true)], summary);
  }

  #[test]
  fn test_iteration_event_serializes_non_finite_values_as_strings() {
    let event: IterationEvent = serde_json::from_value(json!({
      "iteration": 3,
      "n_diff": 0,
      "n_resolved": 0,
      "max_time_change": 0.5,
      "rms_time_change": null,
      "log_lh_seq": null,
      "log_lh_pos": -10.5,
      "log_lh_coal": "inf",
      "log_lh_total": "inf",
      "clock_rate": 0.001,
      "r_squared": null,
    }))
    .unwrap();
    assert_eq!(
      json!({
        "iteration": 3,
        "n_diff": 0,
        "n_resolved": 0,
        "max_time_change": 0.5,
        "rms_time_change": null,
        "log_lh_seq": null,
        "log_lh_pos": -10.5,
        "log_lh_coal": "inf",
        "log_lh_total": "inf",
        "clock_rate": 0.001,
        "r_squared": null,
      }),
      serde_json::to_value(&event).unwrap()
    );
  }

  mod helpers {
    use serde_json::{Value, json};
    use std::fs;
    use std::path::Path;

    pub(super) fn timetree_config(outdir: &Path) -> Value {
      let zika = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/zika/20");
      json!({
        "tree": zika.join("tree.nwk"),
        "metadata": zika.join("metadata.tsv"),
        "alignment": [zika.join("aln.fasta.xz")],
        "max_iter": 2,
        "seed": 7,
        "output_all": outdir,
        "output_selection": ["Tracelog"],
      })
    }

    pub(super) fn tracelog_rows(path: &Path) -> Vec<Vec<String>> {
      let text = fs::read_to_string(path).unwrap();
      let mut lines = text.lines();
      assert_eq!(
        Some("n_diff,n_resolved,max_time_change,rms_time_change,log_lh_seq,log_lh_pos,log_lh_coal,log_lh_total"),
        lines.next(),
        "the tracelog columns are unchanged"
      );
      lines
        .map(|line| {
          line
            .split(',')
            .map(|cell| {
              if cell.is_empty() {
                String::new()
              } else {
                cell.parse::<f64>().unwrap().to_string()
              }
            })
            .collect()
        })
        .collect()
    }
  }
}
