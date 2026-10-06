use crate::port::{PortHeader, PortReply, PortRequest};
use app_commands::app_paths::{AppFolderEnv, AppFolderPath, AppPaths};
use app_commands::app_settings::store::AppSettingsStore;
use app_commands::app_settings::workspace::{OpenedRuns, open_runs_folder};
use app_commands::bridge::error::{ErrorCode, ErrorResponse};
use app_commands::bridge::service::{AppService, Unconfined};
use app_commands::examples_download::{ExampleDownloads, examples_url};
use app_commands::run_limits::RunLimits;
use app_server::error::error_response;
use app_server::routes::{LocalRouters, local_api_routers};
use app_server::state::{DEFAULT_MAX_UPLOAD_SIZE, LocalSettings, ServerConfig};
use axum::Router;
use axum::body::to_bytes;
use axum::response::Response;
use eyre::{Report, WrapErr};
use napi::bindgen_prelude::Uint8Array;
use std::fs;
use std::path::Path;
use std::str::FromStr;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use strum_macros::EnumString;
use tokio::runtime::{Builder, Runtime};
use tokio::task::AbortHandle;
use tokio_stream::StreamExt as _;
use tokio_util::sync::CancellationToken;
use tower::ServiceExt as _;
use treetime_utils::make_internal_report;

const RUNTIME_THREAD_NAME: &str = "treetime-desktop";

pub struct DesktopService {
  routers: LocalRouters,
  runtime: Runtime,
}

impl DesktopService {
  pub fn open(root: &Path, env: &AppFolderEnv) -> Result<Self, Report> {
    let store = Arc::new(AppSettingsStore::open(root)?);
    let settings = store.read()?;
    let mut paths = AppPaths::resolve(root, env, &settings.paths);
    let named_in = settings.paths.runs.is_some().then(|| store.path());
    let OpenedRuns {
      runs,
      error: runs_error,
    } = open_runs_folder(&mut paths, named_in)?;
    fs::create_dir_all(&paths.examples.path).wrap_err_with(|| {
      format!(
        "When creating the examples folder '{}'{}",
        paths.examples.path.display(),
        folder_source(&paths.examples, settings.paths.examples.is_some(), store.path())
      )
    })?;
    let examples = paths.examples.path.clone();
    let app = Arc::new(AppService::new(
      runs,
      examples.clone(),
      Arc::new(Unconfined),
      RunLimits::Settings(Arc::clone(&store)),
    ));
    let routers = local_api_routers(
      app,
      ServerConfig {
        examples_dir: examples.clone(),
        runs_dir: paths.runs.path.clone(),
        max_upload_size: DEFAULT_MAX_UPLOAD_SIZE,
        max_grid_points: None,
        shutdown: CancellationToken::new(),
        settings: Some(LocalSettings {
          store,
          paths,
          runs_error,
          examples: Arc::new(ExampleDownloads::new(examples_url(), &examples)),
        }),
      },
    )?;
    let runtime = Builder::new_multi_thread()
      .enable_all()
      .thread_name(RUNTIME_THREAD_NAME)
      .build()
      .wrap_err("When starting the runtime of the desktop back end")?;
    Ok(Self { routers, runtime })
  }

  pub fn fetch(
    &self,
    request: PortRequest,
    scope: PortScope,
    send: impl Fn(PortReply) -> bool + Send + Sync + 'static,
  ) -> AbortHandle {
    let seq = request.seq;
    let send = Arc::new(send);
    let head_sent = Arc::new(AtomicBool::new(false));
    let response = request
      .into_http()
      .map(|request| self.router(scope).clone().oneshot(request));
    let exchange = {
      let send = Arc::clone(&send);
      let head_sent = Arc::clone(&head_sent);
      self.runtime.spawn(async move {
        let response = match response {
          Ok(response) => match response.await {
            Ok(response) => response,
            Err(infallible) => match infallible {},
          },
          Err(report) => error_response(ErrorResponse::from_report(&report)),
        };
        stream_response(seq, response, &head_sent, send.as_ref()).await;
      })
    };
    let abort = exchange.abort_handle();
    self.runtime.spawn(async move {
      if let Err(err) = exchange.await
        && err.is_panic()
      {
        let error = ErrorResponse::from_panic(err.into_panic().as_ref());
        if head_sent.load(Ordering::SeqCst) {
          send(PortReply::Reset {
            seq,
            message: [error.message, error.causes.join(": ")].join(": "),
          });
        } else {
          stream_response(seq, error_response(error), &head_sent, send.as_ref()).await;
        }
      }
    });
    abort
  }

  pub fn reject_request(&self, seq: u32, message: String) -> Result<Vec<PortReply>, Report> {
    let response = error_response(ErrorResponse {
      code: ErrorCode::InvalidRequest,
      message,
      causes: vec![],
    });
    let (head, body) = response.into_parts();
    let body = self
      .runtime
      .block_on(to_bytes(body, usize::MAX))
      .wrap_err("When reading the body of an error response")?;
    Ok(vec![
      PortReply::Head {
        seq,
        status: head.status.as_u16(),
        headers: PortHeader::list(&head.headers),
      },
      PortReply::Chunk {
        seq,
        data: Uint8Array::new(Vec::from(body)),
      },
      PortReply::End { seq },
    ])
  }

  fn router(&self, scope: PortScope) -> &Router {
    match scope {
      PortScope::Host => &self.routers.host,
      PortScope::Renderer => &self.routers.renderer,
    }
  }
}

/// The router a message port reaches: the host port of the main process also saves run files to chosen paths.
#[derive(Clone, Copy, Debug, PartialEq, Eq, EnumString)]
#[strum(serialize_all = "lowercase")]
pub enum PortScope {
  Host,
  Renderer,
}

impl PortScope {
  pub fn parse(scope: &str) -> Result<Self, Report> {
    Self::from_str(scope).map_err(|err| make_internal_report!("unknown port scope '{scope}': {err}"))
  }
}

pub(crate) async fn stream_response(
  seq: u32,
  response: Response,
  head_sent: &AtomicBool,
  send: &(impl Fn(PortReply) -> bool + Sync),
) {
  let (head, body) = response.into_parts();
  let head = PortReply::Head {
    seq,
    status: head.status.as_u16(),
    headers: PortHeader::list(&head.headers),
  };
  head_sent.store(true, Ordering::SeqCst);
  if !send(head) {
    return;
  }
  let mut frames = body.into_data_stream();
  while let Some(frame) = frames.next().await {
    let reply = match frame {
      Ok(bytes) => PortReply::Chunk {
        seq,
        data: Uint8Array::new(Vec::from(bytes)),
      },
      Err(err) => PortReply::Reset {
        seq,
        message: format!("When reading the response body: {err}"),
      },
    };
    let reset = matches!(reply, PortReply::Reset { .. });
    if !send(reply) || reset {
      return;
    }
  }
  send(PortReply::End { seq });
}

fn folder_source(folder: &AppFolderPath, in_settings: bool, settings_file: &Path) -> String {
  match folder.fixed_by {
    Some(variable) => format!(", set by {variable}"),
    None if in_settings => format!(", named in '{}'", settings_file.display()),
    None => String::new(),
  }
}
