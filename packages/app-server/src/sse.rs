use crate::error::AppError;
use crate::state::AppState;
use app_commands::command::AppCommand;
use app_commands::job::{CancelToken, JobEvent, JobId, JobProgress, JobStarted, TerminalEvent, run_job};
use axum::response::sse::{Event, Sse};
use axum::response::{IntoResponse, Response};
use log::{error, info};
use serde_json::Value;
use std::convert::Infallible;
use std::sync::Arc;
use tokio::sync::mpsc;
use tokio_stream::StreamExt as _;
use tokio_stream::wrappers::UnboundedReceiverStream;
use treetime::cancel::Cancel;

#[allow(
  clippy::disallowed_methods,
  reason = "the synchronous ProgressSink cannot await bounded sends, so SSE buffering does not apply backpressure"
)]
#[expect(
  tail_expr_drop_order,
  reason = "no value in the tail expression depends on drop order"
)]
pub(crate) fn run_command_sse(state: &Arc<AppState>, command: AppCommand, config: Value) -> Response {
  let job_id = JobId::random();
  let handle = match state.jobs.register(job_id.clone()) {
    Ok(handle) => handle,
    Err(err) => return AppError::from(err).into_response(),
  };
  let output_dir = state.config.out_dir.join(job_id.as_str());
  let (tx, rx) = mpsc::unbounded_channel::<JobEvent>();
  drop(tx.send(JobEvent::Started(JobStarted {
    job_id: job_id.clone(),
    command,
  })));

  let job_state = Arc::clone(state);
  let computation = tokio::task::spawn_blocking(move || {
    let cancel = StreamCancel {
      token: handle.token(),
      tx: tx.clone(),
    };
    let progress = JobProgress::new(move |event| drop(tx.send(event)));
    let confine = |config: &mut Value| job_state.paths.confine(command, config, &output_dir);
    run_job(handle.job_id(), command, &config, &confine, &cancel, &progress)
  });

  let stream = async_stream::stream! {
    let mut events = UnboundedReceiverStream::new(rx);
    while let Some(event) = events.next().await {
      yield Ok::<_, Infallible>(sse_event(&event));
    }
    let terminal = match computation.await {
      Ok(terminal) => terminal,
      Err(err) => TerminalEvent::Error {
        job_id: job_id.clone(),
        message: format!("internal error: the computation task failed: {err}"),
        causes: vec![],
      },
    };
    match &terminal {
      TerminalEvent::Ok { .. } => info!("Job {} finished", job_id.as_str()),
      TerminalEvent::Cancelled { .. } => info!("Job {} cancelled", job_id.as_str()),
      TerminalEvent::Error { message, .. } => error!("Job {} failed: {message}", job_id.as_str()),
    }
    yield Ok::<_, Infallible>(sse_event(&JobEvent::Terminal(terminal)));
  };

  Sse::new(stream).into_response()
}

fn sse_event(event: &JobEvent) -> Event {
  let (name, event) = match event {
    JobEvent::Started(data) => ("started", Event::default().json_data(data)),
    JobEvent::Progress(data) => ("progress", Event::default().json_data(data)),
    JobEvent::Log(data) => ("log", Event::default().json_data(data)),
    JobEvent::Iteration(data) => ("iteration", Event::default().json_data(data)),
    JobEvent::Terminal(data) => ("terminal", Event::default().json_data(data)),
  };
  match event {
    Ok(event) => event.event(name),
    Err(err) => Event::default().comment(format!("When serializing a {name} event: {err}")),
  }
}

struct StreamCancel<'a> {
  token: &'a CancelToken,
  tx: mpsc::UnboundedSender<JobEvent>,
}

impl Cancel for StreamCancel<'_> {
  fn is_cancelled(&self) -> bool {
    self.token.is_cancelled() || self.tx.is_closed()
  }
}
