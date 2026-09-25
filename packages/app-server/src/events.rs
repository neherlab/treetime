use crate::error::AppError;
use crate::state::AppState;
use app_commands::job::{JobEvent, JobId};
use app_commands::runs::events::RunEvent;
use axum::response::sse::{Event, KeepAlive, Sse};
use axum::response::{IntoResponse, Response};
use std::convert::Infallible;
use tokio::sync::mpsc;
use tokio_stream::StreamExt as _;
use tokio_stream::wrappers::UnboundedReceiverStream;

#[allow(
  clippy::disallowed_methods,
  reason = "the synchronous event log cannot await bounded sends, so SSE buffering does not apply backpressure"
)]
pub(crate) fn run_events_sse(state: &AppState, id: &JobId, from: usize) -> Result<Response, AppError> {
  let (tx, rx) = mpsc::unbounded_channel::<RunEvent>();
  state.runs.subscribe(
    id,
    from,
    Box::new(move |event: &RunEvent| tx.send(event.clone()).is_ok()),
  )?;
  let stream = UnboundedReceiverStream::new(rx).map(|event| Ok::<_, Infallible>(sse_event(&event)));
  Ok(Sse::new(stream).keep_alive(KeepAlive::default()).into_response())
}

fn sse_event(event: &RunEvent) -> Event {
  let name = match &event.event {
    JobEvent::Started(_) => "started",
    JobEvent::Progress(_) => "progress",
    JobEvent::Log(_) => "log",
    JobEvent::Iteration(_) => "iteration",
    JobEvent::Terminal(_) => "terminal",
  };
  match Event::default().json_data(event) {
    Ok(sse) => sse.event(name).id(event.seq.to_string()),
    Err(err) => Event::default().comment(format!("When serializing a {name} event: {err}")),
  }
}
