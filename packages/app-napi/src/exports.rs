use crate::backend::DesktopService;
use crate::guard::{guarded, to_napi};
use crate::port::{PortReply, PortRequest};
use crate::subscription::EventForwarder;
use app_commands::bridge::operations::OperationRequest;
use app_commands::job::JobId;
use app_commands::runs::errors::parse_request_text;
use eyre::Report;
use napi::bindgen_prelude::{AsyncTask, ToNapiValue, TypeName};
use napi::threadsafe_function::{ThreadsafeFunction, ThreadsafeFunctionCallMode};
use napi::{Env, Status, Task};
use napi_derive::napi;
use std::path::Path;
use std::sync::Arc;
use tokio::task::AbortHandle;
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
  pub fn call(&self, request_json: String) -> AsyncTask<BlockingTask<String>> {
    let service = Arc::clone(&self.service);
    BlockingTask::spawn(move || {
      let request: OperationRequest = parse_request_text(&request_json)?;
      request.handle(service.app())
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
        .app()
        .runs()
        .subscribe(&JobId::parse(&id)?, usize::try_from(from)?, forwarder.subscriber())
    })
    .map_err(|err| to_napi(&err))?;
    Ok(Subscription { forwarder })
  }

  #[napi]
  pub fn fetch(
    &self,
    request: PortRequest,
    on_reply: ThreadsafeFunction<PortReply, (), PortReply, Status, false>,
  ) -> PortExchange {
    let abort = self.service.fetch(request, move |reply| {
      on_reply.call(reply, ThreadsafeFunctionCallMode::NonBlocking) == Status::Ok
    });
    PortExchange { abort }
  }

  #[napi(ts_return_type = "Promise<void>")]
  pub fn save_run_file(&self, request: SaveRunFileRequest) -> AsyncTask<BlockingTask<()>> {
    let service = Arc::clone(&self.service);
    BlockingTask::spawn(move || {
      let SaveRunFileRequest { id, path, destination } = request;
      service.save_run_file(&JobId::parse(&id)?, &path, Path::new(&destination))
    })
  }

  #[napi(ts_return_type = "Promise<void>")]
  pub fn save_run_archive(&self, request: SaveRunArchiveRequest) -> AsyncTask<BlockingTask<()>> {
    let service = Arc::clone(&self.service);
    BlockingTask::spawn(move || {
      let SaveRunArchiveRequest { id, destination } = request;
      service.save_run_archive(&JobId::parse(&id)?, Path::new(&destination))
    })
  }
}

#[napi(object)]
pub struct SaveRunFileRequest {
  pub id: String,
  pub path: String,
  pub destination: String,
}

#[napi(object)]
pub struct SaveRunArchiveRequest {
  pub id: String,
  pub destination: String,
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

type Job<T> = Box<dyn FnOnce() -> Result<T, Report> + Send>;

pub struct BlockingTask<T> {
  job: Option<Job<T>>,
}

impl<T: ToNapiValue + TypeName + Send + 'static> BlockingTask<T> {
  fn spawn(job: impl FnOnce() -> Result<T, Report> + Send + 'static) -> AsyncTask<Self> {
    AsyncTask::new(Self {
      job: Some(Box::new(job)),
    })
  }
}

impl<T: ToNapiValue + TypeName + Send + 'static> Task for BlockingTask<T> {
  type Output = T;
  type JsValue = T;

  fn compute(&mut self) -> napi::Result<Self::Output> {
    let job = self.job.take();
    guarded(|| job.map_or_else(|| Err(make_report!("the task has already run")), |job| job()))
      .map_err(|err| to_napi(&err))
  }

  fn resolve(&mut self, _env: Env, output: T) -> napi::Result<T> {
    Ok(output)
  }
}
