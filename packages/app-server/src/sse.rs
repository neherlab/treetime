use crate::error::AppError;
use axum::response::sse::{Event, Sse};
use axum::response::{IntoResponse, Response};
use eyre::Report;
use log::{error, info};
use serde::Serialize;
use serde_json::Value;
use std::convert::Infallible;
use std::path::Path;
use tokio::sync::mpsc;
use tokio_stream::StreamExt as _;
use tokio_stream::wrappers::UnboundedReceiverStream;
use treetime::cancel::{Cancel, CancelledError};
use treetime::progress::{LogEvent, LogLevel, ProgressSink};
use treetime_schema::ProgressEvent;

enum SinkEvent {
  Progress(ProgressEvent),
  Log(LogEvent),
}

struct ChannelProgress {
  tx: mpsc::UnboundedSender<SinkEvent>,
}

impl ChannelProgress {
  fn new(tx: mpsc::UnboundedSender<SinkEvent>) -> Self {
    Self { tx }
  }
}

impl ProgressSink for ChannelProgress {
  fn report(&self, stage: &str, fraction: f64, message: &str) {
    drop(self.tx.send(SinkEvent::Progress(ProgressEvent {
      stage: stage.to_owned(),
      fraction,
      message: message.to_owned(),
    })));
  }

  fn log(&self, level: LogLevel, message: &str) {
    drop(self.tx.send(SinkEvent::Log(LogEvent {
      level,
      message: message.to_owned(),
    })));
  }

  fn log_enabled(&self, _level: LogLevel) -> bool {
    true
  }
}

impl Cancel for ChannelProgress {
  fn is_cancelled(&self) -> bool {
    self.tx.is_closed()
  }
}

#[allow(
  clippy::disallowed_methods,
  clippy::expect_used,
  reason = "the synchronous ProgressSink cannot await bounded sends, so SSE buffering does not apply backpressure; event serialization is infallible"
)]
#[expect(
  tail_expr_drop_order,
  reason = "no value in the tail expression depends on drop order"
)]
fn sse_response<F>(run_fn: F) -> Response
where
  F: FnOnce(&dyn Cancel, &dyn ProgressSink) -> Result<Value, Report> + Send + 'static,
{
  let (tx, rx) = mpsc::unbounded_channel::<SinkEvent>();

  let computation = tokio::task::spawn_blocking(move || {
    let progress = ChannelProgress::new(tx);
    run_fn(&progress, &progress)
  });

  let stream = async_stream::stream! {
    let mut rx_stream = UnboundedReceiverStream::new(rx);
    while let Some(event) = rx_stream.next().await {
      match event {
        SinkEvent::Progress(p) => {
          yield Ok::<_, Infallible>(
            Event::default()
              .event("progress")
              .json_data(p)
              .expect("ProgressEvent serialization"),
          );
        },
        SinkEvent::Log(l) => {
          yield Ok::<_, Infallible>(
            Event::default()
              .event("log")
              .json_data(l)
              .expect("LogEvent serialization"),
          );
        },
      }
    }

    match computation.await {
      Ok(Ok(value)) => {
        yield Ok::<_, Infallible>(
          Event::default()
            .event("result")
            .json_data(value)
            .expect("result serialization"),
        );
      },
      Ok(Err(err)) if err.downcast_ref::<CancelledError>().is_some() => {
        info!("Computation cancelled by client");
      },
      Ok(Err(err)) => {
        error!("Computation failed: {err:?}");
        yield Ok::<_, Infallible>(
          Event::default()
            .event("result")
            .json_data(serde_json::json!({ "error": format!("{err:#}") }))
            .expect("result serialization"),
        );
      },
      Err(err) => {
        error!("Computation panicked: {err}");
        yield Ok::<_, Infallible>(
          Event::default()
            .event("result")
            .json_data(serde_json::json!({ "error": format!("{err}") }))
            .expect("result serialization"),
        );
      },
    }
  };

  Sse::new(stream).into_response()
}

pub(crate) fn handle_command<S, T>(
  mut body: Value,
  out_dir: &Path,
  command: fn(&S, &dyn Cancel, &dyn ProgressSink) -> Result<T, Report>,
) -> Response
where
  S: serde::de::DeserializeOwned + Send + 'static,
  T: Serialize + Send + 'static,
{
  if let Some(obj) = body.as_object_mut() {
    let client_outdir = obj.get("outdir").and_then(Value::as_str).unwrap_or_default().to_owned();
    let resolved = out_dir.join(client_outdir);
    obj.insert(
      "outdir".to_owned(),
      Value::String(resolved.to_string_lossy().into_owned()),
    );
  }
  let args: S = match serde_json::from_value(body) {
    Ok(args) => args,
    Err(err) => return AppError::from(err).into_response(),
  };
  sse_response(move |cancel, progress| {
    let result = command(&args, cancel, progress)?;
    serde_json::to_value(result).map_err(Report::from)
  })
}
