use crate::backend::DesktopService;
use crate::guard::{guarded, to_napi};
use crate::subscription::EventForwarder;
use app_commands::bridge::operations::DesktopRequest;
use app_commands::job::JobId;
use eyre::Report;
use napi::bindgen_prelude::AsyncTask;
use napi::threadsafe_function::{ThreadsafeFunction, ThreadsafeFunctionCallMode};
use napi::{Env, Status, Task};
use napi_derive::napi;
use std::path::Path;
use std::sync::Arc;
use treetime_utils::make_report;

type EventCallback = Arc<ThreadsafeFunction<String, ()>>;

type EventSink = Box<dyn Fn(String) -> bool + Send>;

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
  pub fn new(runs_dir: String) -> napi::Result<Self> {
    let service = guarded(|| DesktopService::open(Path::new(&runs_dir))).map_err(|err| to_napi(&err))?;
    Ok(Self {
      service: Arc::new(service),
    })
  }

  #[napi(ts_return_type = "Promise<string>")]
  pub fn call(&self, request_json: String) -> AsyncTask<JsonTask> {
    let service = Arc::clone(&self.service);
    JsonTask::spawn(move || {
      let request: DesktopRequest = serde_json::from_str(&request_json)?;
      request.handle(service.as_ref())
    })
  }

  #[napi(ts_args_type = "id: string, from: number, onEvent: (err: Error | null, eventJson: string) => void")]
  pub fn subscribe(&self, id: String, from: u32, on_event: EventCallback) -> napi::Result<Subscription> {
    let send: EventSink =
      Box::new(move |json: String| on_event.call(Ok(json), ThreadsafeFunctionCallMode::NonBlocking) == Status::Ok);
    let forwarder = EventForwarder::new(send);
    guarded(|| {
      self
        .service
        .runs()
        .subscribe(&JobId::parse(&id)?, usize::try_from(from)?, forwarder.subscriber())
    })
    .map_err(|err| to_napi(&err))?;
    Ok(Subscription { forwarder })
  }

  #[napi(ts_return_type = "Promise<string>")]
  pub fn save_run_file(&self, id: String, path: String, destination: String) -> AsyncTask<JsonTask> {
    let service = Arc::clone(&self.service);
    JsonTask::spawn(move || {
      service.save_run_file(&JobId::parse(&id)?, &path, Path::new(&destination))?;
      Ok(serde_json::to_string(&())?)
    })
  }

  #[napi(ts_return_type = "Promise<string>")]
  pub fn save_run_archive(&self, id: String, destination: String) -> AsyncTask<JsonTask> {
    let service = Arc::clone(&self.service);
    JsonTask::spawn(move || {
      service.save_run_archive(&JobId::parse(&id)?, Path::new(&destination))?;
      Ok(serde_json::to_string(&())?)
    })
  }
}

#[napi]
pub struct Subscription {
  forwarder: EventForwarder<EventSink>,
}

#[napi]
impl Subscription {
  #[napi]
  pub fn unsubscribe(&self) {
    self.forwarder.close();
  }
}

type Job = Box<dyn FnOnce() -> Result<String, Report> + Send>;

pub struct JsonTask {
  job: Option<Job>,
}

impl JsonTask {
  fn spawn(job: impl FnOnce() -> Result<String, Report> + Send + 'static) -> AsyncTask<Self> {
    AsyncTask::new(Self {
      job: Some(Box::new(job)),
    })
  }
}

impl Task for JsonTask {
  type Output = String;
  type JsValue = String;

  fn compute(&mut self) -> napi::Result<Self::Output> {
    let job = self.job.take();
    guarded(|| job.map_or_else(|| Err(make_report!("the task has already run")), |job| job()))
      .map_err(|err| to_napi(&err))
  }

  fn resolve(&mut self, _env: Env, output: String) -> napi::Result<String> {
    Ok(output)
  }
}
