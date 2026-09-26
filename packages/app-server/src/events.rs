use crate::api::response::TypedSse;
use crate::error::AppError;
use crate::state::AppState;
use app_commands::job::JobId;
use app_commands::runs::events::RunEvent;
use axum::response::sse::Event;
use tokio::sync::mpsc;
use tokio_stream::wrappers::UnboundedReceiverStream;

#[allow(
  clippy::disallowed_methods,
  reason = "the synchronous event log cannot await bounded sends, so SSE buffering does not apply backpressure"
)]
pub(crate) fn run_events_sse(state: &AppState, id: &JobId, from: usize) -> Result<TypedSse<RunEvent>, AppError> {
  let (tx, rx) = mpsc::unbounded_channel::<RunEvent>();
  state.runs.subscribe(
    id,
    from,
    Box::new(move |event: &RunEvent| tx.send(event.clone()).is_ok()),
  )?;
  Ok(TypedSse::new(UnboundedReceiverStream::new(rx), sse_event))
}

fn sse_event(event: &RunEvent) -> Event {
  let name: &'static str = (&event.event).into();
  match Event::default().json_data(event) {
    Ok(sse) => sse.event(name).id(event.seq.to_string()),
    Err(err) => Event::default().comment(format!("When serializing a {name} event: {err}")),
  }
}
