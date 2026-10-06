use crate::api::response::TypedSse;
use crate::error::AppError;
use crate::state::AppState;
use app_commands::job::JobId;
use app_commands::runs::app_events::AppEvent;
use app_commands::runs::events::RunEvent;
use axum::response::sse::Event;
use deser::Serialize;
use std::convert;
use tokio::sync::mpsc;
use tokio_stream::wrappers::UnboundedReceiverStream;
use tokio_stream::{self as stream, Stream, StreamExt as _};
use tokio_util::sync::CancellationToken;
use treetime_utils::io::json::{JsonPretty, json_write_str};

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
  let events = until_shutdown(UnboundedReceiverStream::new(rx), state.config.shutdown.clone());
  Ok(TypedSse::new(events, run_sse_event))
}

#[allow(
  clippy::disallowed_methods,
  reason = "the synchronous event log cannot await bounded sends, so SSE buffering does not apply backpressure"
)]
pub(crate) fn app_events_sse(state: &AppState, from: Option<usize>) -> TypedSse<AppEvent> {
  let (tx, rx) = mpsc::unbounded_channel::<AppEvent>();
  state
    .runs
    .app_events()
    .subscribe(from, Box::new(move |event: &AppEvent| tx.send(event.clone()).is_ok()));
  let events = until_shutdown(UnboundedReceiverStream::new(rx), state.config.shutdown.clone());
  TypedSse::new(events, app_sse_event)
}

fn until_shutdown<T: Send + 'static>(
  items: impl Stream<Item = T> + Send + 'static,
  shutdown: CancellationToken,
) -> impl Stream<Item = T> + Send + 'static {
  let stop = stream::once(())
    .then(move |()| shutdown.clone().cancelled_owned())
    .map(|()| None);
  items
    .map(Some)
    .chain(stream::once(None))
    .merge(stop)
    .map_while(convert::identity)
}

fn run_sse_event(event: &RunEvent) -> Event {
  sse_event((&event.event).into(), event.seq, event)
}

fn app_sse_event(event: &AppEvent) -> Event {
  sse_event((&event.change).into(), event.seq, event)
}

fn sse_event(name: &'static str, seq: usize, event: &impl Serialize) -> Event {
  match json_write_str(event, JsonPretty(false)) {
    Ok(data) => Event::default().data(data).event(name).id(seq.to_string()),
    Err(err) => Event::default().comment(format!("When serializing a {name} event: {err}")),
  }
}
