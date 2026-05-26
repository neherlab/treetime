use crate::error::AppError;
use app_api::progress::ProgressSink;
use axum::response::sse::{Event, Sse};
use axum::response::{IntoResponse, Response};
use eyre::Report;
use serde_json::Value;
use std::convert::Infallible;
use tokio::sync::mpsc;
use tokio_stream::StreamExt as _;
use tokio_stream::wrappers::UnboundedReceiverStream;
use treetime::progress::{LogLevel, ProgressEvent};

pub struct ChannelProgress {
  tx: mpsc::UnboundedSender<ProgressEvent>,
}

impl ChannelProgress {
  fn new(tx: mpsc::UnboundedSender<ProgressEvent>) -> Self {
    Self { tx }
  }
}

impl ProgressSink for ChannelProgress {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    drop(self.tx.send(ProgressEvent {
      stage: stage.to_owned(),
      fraction,
      message: message.to_owned(),
    }));
  }

  fn log(&self, _level: LogLevel, _message: &str) {}
}

#[allow(tail_expr_drop_order)]
pub fn sse_response<F>(run_fn: F) -> Response
where
  F: FnOnce(&dyn ProgressSink) -> Result<Value, Report> + Send + 'static,
{
  let (tx, rx) = mpsc::unbounded_channel::<ProgressEvent>();

  let computation = tokio::task::spawn_blocking(move || {
    let progress = ChannelProgress::new(tx);
    run_fn(&progress)
  });

  let stream = async_stream::stream! {
    let mut rx_stream = UnboundedReceiverStream::new(rx);
    while let Some(event) = rx_stream.next().await {
      yield Ok::<_, Infallible>(
        Event::default()
          .event("progress")
          .json_data(event)
          .expect("ProgressEvent serialization"),
      );
    }

    let result_value = match computation.await {
      Ok(Ok(value)) => value,
      Ok(Err(err)) => serde_json::json!({ "error": format!("{err:?}") }),
      Err(err) => serde_json::json!({ "error": format!("computation panicked: {err}") }),
    };
    yield Ok::<_, Infallible>(
      Event::default()
        .event("result")
        .json_data(result_value)
        .expect("result serialization"),
    );
  };

  Sse::new(stream).into_response()
}

#[allow(clippy::result_large_err)]
pub fn parse_args<A: serde::de::DeserializeOwned>(body: Value) -> Result<A, Response> {
  serde_json::from_value(body).map_err(|err| AppError::from(err).into_response())
}
