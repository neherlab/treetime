#[cfg(test)]
mod tests {
  use crate::backend::web_router;
  use crate::exchange::PortExchange;
  use crate::jobs::{JobExecutor, JobOrder};
  use app_backend::runs::events::read_events;
  use app_backend::runs::lifecycle::{JobInput, RunJob};
  use axum::Router;
  use helpers::{call, data_file, started_run};
  use pretty_assertions::assert_eq;
  use serde_json::json;
  use std::fs;
  use std::path::PathBuf;
  use treetime_utils::assert_error;
  use treetime_utils::io::json::json_read_str;

  #[tokio::test]
  async fn test_jobs_start_order_names_each_input_file_once() {
    let run = started_run().await;
    let orders = run.jobs.take_orders().unwrap();
    let [order] = orders.as_slice() else {
      panic!("expected one order, got {orders:?}");
    };
    let job: RunJob = json_read_str(order.job.as_deref().unwrap()).unwrap();
    let paths = vec![data_file("aln.fasta.xz"), data_file("tree.nwk")];
    assert_eq!(
      (
        "start",
        run.id.as_str(),
        paths.iter().map(|path| path.display().to_string()).collect::<Vec<_>>(),
        paths,
      ),
      (
        order.kind.as_str(),
        order.run_id.as_str(),
        order.input_paths.clone(),
        job.inputs.iter().map(|input: &JobInput| input.path.clone()).collect(),
      )
    );
    assert_eq!(Vec::<JobOrder>::new(), run.jobs.take_orders().unwrap());
  }

  #[tokio::test]
  async fn test_jobs_cancel_ends_the_run_cancelled_and_refuses_a_later_completion() {
    let run = started_run().await;
    run.jobs.take_orders().unwrap();
    let (cancel_status, cancel_body) = call(&run.router, "POST", &format!("/api/runs/{}/cancel", run.id)).await;
    let orders = run.jobs.take_orders().unwrap();
    let failed = run.jobs.fail(&run.id, vec![]).unwrap();
    let outcome = json!({ "terminal": { "status": "cancelled", "job_id": run.id }, "warnings": [] }).to_string();
    let completed = run.jobs.complete(&run.id, &outcome).unwrap();
    let (_, record) = call(&run.router, "GET", &format!("/api/runs/{}", run.id)).await;
    assert_eq!(
      (
        200,
        json!({ "cancelled": true }),
        vec![("cancel", run.id.clone())],
        true,
        false,
        json!("cancelled")
      ),
      (
        cancel_status,
        cancel_body,
        orders
          .iter()
          .map(|order| (order.kind.as_str(), order.run_id.clone()))
          .collect::<Vec<_>>(),
        failed,
        completed,
        record["status"].clone()
      )
    );
  }

  #[tokio::test]
  async fn test_jobs_failure_records_the_last_line_as_message_and_the_earlier_lines_as_causes() {
    let run = started_run().await;
    let failed = run
      .jobs
      .fail(
        &run.id,
        vec![
          "thread 'main' panicked at src/lib.rs:1:1: index out of bounds".to_owned(),
          "unreachable".to_owned(),
        ],
      )
      .unwrap();
    let (_, record) = call(&run.router, "GET", &format!("/api/runs/{}", run.id)).await;
    assert_eq!(
      (
        true,
        json!("error"),
        json!({
          "message": "unreachable",
          "causes": ["thread 'main' panicked at src/lib.rs:1:1: index out of bounds"],
        })
      ),
      (failed, record["status"].clone(), record["error"].clone())
    );
  }

  #[tokio::test]
  async fn test_jobs_failure_without_lines_says_the_worker_stopped_without_a_message() {
    let run = started_run().await;
    run.jobs.fail(&run.id, vec![]).unwrap();
    let (_, record) = call(&run.router, "GET", &format!("/api/runs/{}", run.id)).await;
    assert_eq!(
      json!({ "message": "the run worker stopped without an error message", "causes": [] }),
      record["error"]
    );
  }

  #[tokio::test]
  async fn test_jobs_output_paths_must_stay_inside_the_output_folder() {
    let run = started_run().await;
    assert_error!(
      run.jobs.begin_output(&run.id, "../run.json"),
      "file path `../run.json` must name a file inside the run's output folder"
    );
    assert_error!(
      run.jobs.begin_output(&run.id, "/treetime/runs/x"),
      "file path `/treetime/runs/x` must name a file inside the run's output folder"
    );
  }

  #[tokio::test]
  async fn test_jobs_an_empty_output_file_exists_after_its_begin() {
    let run = started_run().await;
    run.jobs.begin_output(&run.id, "nested/empty.txt").unwrap();
    run.jobs.begin_output(&run.id, "nested/tree.nwk").unwrap();
    run.jobs.append_output(&run.id, "nested/tree.nwk", b"(A,").unwrap();
    run.jobs.append_output(&run.id, "nested/tree.nwk", b"B);").unwrap();
    let out_dir = run.runs_dir.path().join(&run.id).join("out").join("nested");
    assert_eq!(
      (String::new(), "(A,B);".to_owned()),
      (
        fs::read_to_string(out_dir.join("empty.txt")).unwrap(),
        fs::read_to_string(out_dir.join("tree.nwk")).unwrap()
      )
    );
  }

  #[tokio::test]
  async fn test_jobs_a_worker_cannot_report_the_start_or_the_end_of_a_run() {
    let run = started_run().await;
    let started = json!({ "type": "started", "data": { "job_id": run.id, "command": "ancestral" } }).to_string();
    assert_error!(
      run.jobs.append_event(&run.id, &started),
      format!(
        "a run worker cannot report a `started` event of run `{}`; only the back end starts and ends a run",
        run.id
      )
    );
  }

  #[tokio::test]
  async fn test_jobs_progress_events_reach_the_event_log_of_the_run() {
    let run = started_run().await;
    let progress = json!({ "type": "log", "data": { "level": "info", "message": "half way" } }).to_string();
    run.jobs.append_event(&run.id, &progress).unwrap();
    let events = read_events(&run.runs_dir.path().join(&run.id).join("events.jsonl"), 0).unwrap();
    assert_eq!(
      vec!["started", "log"],
      events
        .iter()
        .map(|event| <&str>::from(&event.event))
        .collect::<Vec<_>>()
    );
  }

  #[tokio::test]
  async fn test_jobs_a_worker_reads_only_the_inputs_of_its_run() {
    let run = started_run().await;
    assert_error!(
      run.jobs.read_input(&run.id, "/etc/passwd"),
      format!("'/etc/passwd' is not an input file of run `{}`", run.id)
    );
  }

  #[tokio::test]
  async fn test_jobs_a_run_that_is_not_waiting_refuses_worker_messages() {
    let jobs = JobExecutor::default();
    assert_error!(
      jobs.begin_output("r1", "tree.nwk"),
      "run `r1` is not computing in a run worker"
    );
  }

  mod helpers {
    use super::{JobExecutor, PathBuf, PortExchange, Router, web_router};
    use app_backend::bridge::executor::RunExecutor;
    use app_server::port::port_request;
    use app_server::state::{DEFAULT_MAX_UPLOAD_SIZE, ServerConfig};
    use serde_json::{Value, json};
    use std::sync::Arc;
    use tempfile::{TempDir, tempdir};
    use tokio_util::sync::CancellationToken;

    pub(super) struct StartedRun {
      pub router: Router,
      pub jobs: Arc<JobExecutor>,
      pub id: String,
      pub runs_dir: TempDir,
    }

    pub(super) async fn started_run() -> StartedRun {
      let runs_dir = tempdir().unwrap();
      let jobs = Arc::new(JobExecutor::default());
      let executor: Arc<dyn RunExecutor> = Arc::<JobExecutor>::clone(&jobs);
      let router = web_router(
        ServerConfig {
          examples_dir: data_dir(),
          runs_dir: runs_dir.path().to_owned(),
          max_upload_size: DEFAULT_MAX_UPLOAD_SIZE,
          max_grid_points: None,
          shutdown: CancellationToken::new(),
          settings: None,
        },
        executor,
      )
      .unwrap();
      let alignment = data_file("aln.fasta.xz");
      let body = json!({
        "command": "ancestral",
        "config": { "tree": data_file("tree.nwk"), "alignment": [alignment, alignment] },
      });
      let (status, record) = send(&router, "POST", "/api/runs", Some(body)).await;
      assert_eq!(200, status, "{record}");
      StartedRun {
        router,
        jobs,
        id: record["id"].as_str().unwrap().to_owned(),
        runs_dir,
      }
    }

    pub(super) async fn call(router: &Router, method: &str, url: &str) -> (u16, Value) {
      send(router, method, url, None).await
    }

    pub(super) fn data_file(name: &str) -> PathBuf {
      data_dir().join("flu/h3n2/20").join(name)
    }

    fn data_dir() -> PathBuf {
      PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../data")
        .canonicalize()
        .unwrap()
    }

    async fn send(router: &Router, method: &str, url: &str, body: Option<Value>) -> (u16, Value) {
      let headers = vec![("content-type".to_owned(), "application/json".to_owned())];
      let body = body.map(|body| body.to_string().into_bytes());
      let exchange = PortExchange::start(router.clone(), port_request(method, url, headers, body));
      let head = exchange.head().await.unwrap();
      let mut bytes = vec![];
      while let Some(chunk) = exchange.chunk().await.unwrap() {
        bytes.extend_from_slice(&chunk);
      }
      (head.status, serde_json::from_slice(&bytes).unwrap())
    }
  }
}
