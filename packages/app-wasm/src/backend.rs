use crate::exchange::{Exchange, PortExchange};
use app_commands::runs::errors::invalid;
use app_server::port::port_request;
use app_server::routes::api_router;
use app_server::state::{DEFAULT_MAX_UPLOAD_SIZE, ServerConfig, server_service};
use axum::Router;
use eyre::{Report, WrapErr};
use itertools::izip;
use std::fs;
use std::path::PathBuf;
use tokio_util::sync::CancellationToken;
use treetime_utils::error::report_to_string;
use wasm_bindgen::prelude::{JsError, wasm_bindgen};

const EXAMPLES_DIR: &str = "/treetime/examples";

const RUNS_DIR: &str = "/treetime/runs";

#[wasm_bindgen]
pub struct Backend {
  router: Router,
}

#[wasm_bindgen]
impl Backend {
  #[wasm_bindgen(constructor)]
  pub fn new() -> Result<Self, JsError> {
    create_examples_folder()
      .and_then(|()| web_router(browser_config()))
      .map(|router| Self { router })
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
}

pub(crate) fn web_router(config: ServerConfig) -> Result<Router, Report> {
  let (router, _) = api_router(server_service(&config)?, config)?;
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
