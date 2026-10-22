use crate::chunks::ChunkReader;
use crate::exchange::{Exchange, PortExchange};
use crate::jobs::{JobExecutor, JobOrder};
use app_backend::bridge::executor::RunExecutor;
use app_backend::runs::errors::invalid;
use app_server::port::port_request;
use app_server::routes::api_router;
use app_server::state::{DEFAULT_MAX_UPLOAD_SIZE, ServerConfig, server_service};
use axum::Router;
use eyre::{Report, WrapErr};
use itertools::izip;
use std::fs;
use std::path::PathBuf;
use std::sync::Arc;
use tokio_util::sync::CancellationToken;
use treetime_utils::error::report_to_string;
use wasm_bindgen::prelude::{JsError, wasm_bindgen};

const EXAMPLES_DIR: &str = "/treetime/examples";

const RUNS_DIR: &str = "/treetime/runs";

#[wasm_bindgen]
pub struct Backend {
  router: Router,
  pub(crate) jobs: Arc<JobExecutor>,
}

#[wasm_bindgen]
impl Backend {
  #[wasm_bindgen(constructor)]
  pub fn new() -> Result<Self, JsError> {
    let jobs = Arc::new(JobExecutor::default());
    let executor: Arc<dyn RunExecutor> = Arc::<JobExecutor>::clone(&jobs);
    create_examples_folder()
      .and_then(|()| web_router(browser_config(), executor))
      .map(|router| Self { router, jobs })
      .map_err(|report| JsError::new(&report_to_string(&report)))
  }

  pub fn fetch(
    &self,
    method: &str,
    url: &str,
    #[wasm_bindgen(js_name = headerNames)] header_names: Vec<String>,
    #[wasm_bindgen(js_name = headerValues)] header_values: Vec<String>,
    body: Option<Vec<u8>>,
  ) -> Exchange {
    Exchange {
      inner: PortExchange::start(
        self.router.clone(),
        header_list(header_names, header_values).and_then(|headers| port_request(method, url, headers, body)),
      ),
    }
  }

  #[wasm_bindgen(js_name = takeJobOrders)]
  pub fn take_job_orders(&self) -> Result<Vec<JobOrder>, JsError> {
    self.jobs.take_orders().map_err(|report| js_error(&report))
  }

  #[wasm_bindgen(js_name = readRunInput)]
  pub fn read_run_input(
    &self,
    #[wasm_bindgen(js_name = runId)] run_id: &str,
    path: &str,
  ) -> Result<ChunkReader, JsError> {
    self
      .jobs
      .read_input(run_id, path)
      .map(ChunkReader::from)
      .map_err(|report| js_error(&report))
  }

  #[wasm_bindgen(js_name = appendRunEvent)]
  pub fn append_run_event(
    &self,
    #[wasm_bindgen(js_name = runId)] run_id: &str,
    #[wasm_bindgen(js_name = eventJson)] event_json: &str,
  ) -> Result<(), JsError> {
    self
      .jobs
      .append_event(run_id, event_json)
      .map_err(|report| js_error(&report))
  }

  #[wasm_bindgen(js_name = beginRunOutput)]
  pub fn begin_run_output(
    &self,
    #[wasm_bindgen(js_name = runId)] run_id: &str,
    #[wasm_bindgen(js_name = relativePath)] relative_path: &str,
  ) -> Result<(), JsError> {
    self
      .jobs
      .begin_output(run_id, relative_path)
      .map_err(|report| js_error(&report))
  }

  #[wasm_bindgen(js_name = appendRunOutput)]
  pub fn append_run_output(
    &self,
    #[wasm_bindgen(js_name = runId)] run_id: &str,
    #[wasm_bindgen(js_name = relativePath)] relative_path: &str,
    chunk: &[u8],
  ) -> Result<(), JsError> {
    self
      .jobs
      .append_output(run_id, relative_path, chunk)
      .map_err(|report| js_error(&report))
  }

  #[wasm_bindgen(js_name = completeRun)]
  pub fn complete_run(
    &self,
    #[wasm_bindgen(js_name = runId)] run_id: &str,
    #[wasm_bindgen(js_name = outcomeJson)] outcome_json: &str,
  ) -> Result<bool, JsError> {
    self
      .jobs
      .complete(run_id, outcome_json)
      .map_err(|report| js_error(&report))
  }

  #[wasm_bindgen(js_name = failRun)]
  pub fn fail_run(&self, #[wasm_bindgen(js_name = runId)] run_id: &str, lines: Vec<String>) -> Result<bool, JsError> {
    self.jobs.fail(run_id, lines).map_err(|report| js_error(&report))
  }
}

pub(crate) fn web_router(config: ServerConfig, executor: Arc<dyn RunExecutor>) -> Result<Router, Report> {
  let (router, _) = api_router(server_service(&config, executor)?, config)?;
  Ok(router)
}

fn create_examples_folder() -> Result<(), Report> {
  fs::create_dir_all(EXAMPLES_DIR).wrap_err_with(|| format!("When creating the examples folder '{EXAMPLES_DIR}'"))
}

fn browser_config() -> ServerConfig {
  ServerConfig {
    examples_dir: PathBuf::from(EXAMPLES_DIR),
    runs_dir: PathBuf::from(RUNS_DIR),
    max_upload_size: DEFAULT_MAX_UPLOAD_SIZE,
    max_grid_points: None,
    shutdown: CancellationToken::new(),
    settings: None,
  }
}

pub(crate) fn header_list(names: Vec<String>, values: Vec<String>) -> Result<Vec<(String, String)>, Report> {
  if names.len() != values.len() {
    return Err(invalid(format!(
      "a request has {} header names and {} header values",
      names.len(),
      values.len()
    )));
  }
  Ok(izip!(names, values).collect())
}

fn js_error(report: &Report) -> JsError {
  JsError::new(&report_to_string(report))
}
