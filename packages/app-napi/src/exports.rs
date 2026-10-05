use crate::backend::{DesktopService, PortScope};
use crate::guard::{guarded, to_napi};
use crate::port::{PortReply, PortRequest};
use app_commands::app_paths::{AppFolderEnv, AppPaths, app_root};
use app_commands::app_settings::settings::UiTheme;
use app_commands::app_settings::store::AppSettingsStore;
use eyre::Report;
use napi::Status;
use napi::threadsafe_function::{ThreadsafeFunction, ThreadsafeFunctionCallMode};
use napi_derive::napi;
use serde_json::Value;
use std::path::Path;
use std::sync::Arc;
use tokio::task::AbortHandle;
use treetime_utils::{make_internal_report, make_report};

#[napi]
pub struct Backend {
  service: Arc<DesktopService>,
}

#[napi]
#[allow(
  clippy::needless_pass_by_value,
  reason = "napi passes JavaScript values as owned arguments; the napi macro re-emits the item, so expect cannot track it"
)]
impl Backend {
  #[napi(constructor)]
  pub fn new() -> napi::Result<Self> {
    let service =
      guarded(|| DesktopService::open(&app_root()?, &AppFolderEnv::from_env()?)).map_err(|err| to_napi(&err))?;
    Ok(Self {
      service: Arc::new(service),
    })
  }

  #[napi(ts_args_type = "request: PortRequest, scope: 'host' | 'renderer', onReply: ((arg: PortReply) => void)")]
  pub fn fetch(
    &self,
    request: PortRequest,
    scope: String,
    on_reply: ThreadsafeFunction<PortReply, (), PortReply, Status, false>,
  ) -> napi::Result<PortExchange> {
    let scope = guarded(|| PortScope::parse(&scope)).map_err(|err| to_napi(&err))?;
    let abort = self.service.fetch(request, scope, move |reply| {
      on_reply.call(reply, ThreadsafeFunctionCallMode::NonBlocking) == Status::Ok
    });
    Ok(PortExchange { abort })
  }

  #[napi]
  pub fn reject_request(&self, seq: u32, message: String) -> napi::Result<Vec<PortReply>> {
    guarded(|| self.service.reject_request(seq, message)).map_err(|err| to_napi(&err))
  }
}

#[napi]
pub fn app_startup() -> napi::Result<AppStartup> {
  guarded(|| {
    let root = app_root()?;
    let settings = AppSettingsStore::open(&root)?.read()?;
    let paths = AppPaths::resolve(&root, &AppFolderEnv::from_env()?, &settings.paths);
    let Value::String(theme) = serde_json::to_value(settings.ui.theme.unwrap_or(UiTheme::System))? else {
      return Err(make_internal_report!("a theme serializes to a JSON string"));
    };
    Ok(AppStartup {
      profile_dir: path_string(&paths.profile.path)?,
      logs_dir: path_string(&paths.logs.path)?,
      theme,
    })
  })
  .map_err(|err| to_napi(&err))
}

#[napi(object)]
pub struct AppStartup {
  pub profile_dir: String,
  pub logs_dir: String,
  pub theme: String,
}

#[napi]
pub struct PortExchange {
  abort: AbortHandle,
}

#[napi]
impl PortExchange {
  #[napi]
  pub fn abort(&self) {
    self.abort.abort();
  }
}

fn path_string(path: &Path) -> Result<String, Report> {
  path
    .to_str()
    .map(str::to_owned)
    .ok_or_else(|| make_report!("the path '{}' is not valid UTF-8", path.display()))
}
