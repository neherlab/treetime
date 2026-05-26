use crate::error::AppError;
use app_api::progress::ProgressSink;
use axum::response::sse::{Event, Sse};
use axum::response::{IntoResponse, Response};
use eyre::Report;
use log::error;
use serde::Serialize;
use serde_json::Value;
use std::convert::Infallible;
use std::path::Path;
use tokio::sync::mpsc;
use tokio_stream::StreamExt as _;
use tokio_stream::wrappers::UnboundedReceiverStream;
use treetime::progress::{LogEvent, LogLevel, ProgressEvent};

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

#[allow(tail_expr_drop_order)]
fn sse_response<F>(run_fn: F) -> Response
where
  F: FnOnce(&dyn ProgressSink) -> Result<Value, Report> + Send + 'static,
{
  let (tx, rx) = mpsc::unbounded_channel::<SinkEvent>();

  let computation = tokio::task::spawn_blocking(move || {
    let progress = ChannelProgress::new(tx);
    run_fn(&progress)
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

    let result_value = match computation.await {
      Ok(Ok(value)) => value,
      Ok(Err(err)) => {
        error!("Computation failed: {err:?}");
        serde_json::json!({ "error": format!("{err:#}") })
      },
      Err(err) => {
        error!("Computation panicked: {err}");
        serde_json::json!({ "error": format!("{err}") })
      },
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

pub fn handle_command<S, R, T>(
  mut body: Value,
  out_dir: &Path,
  command: fn(&R, &dyn ProgressSink) -> Result<T, Report>,
) -> Response
where
  S: serde::de::DeserializeOwned + Into<R> + Send + 'static,
  R: Send + 'static,
  T: Serialize + Send + 'static,
{
  if let Some(obj) = body.as_object_mut() {
    let client_outdir = obj
      .get("outdir")
      .and_then(Value::as_str)
      .unwrap_or_default()
      .to_owned();
    let resolved = out_dir.join(client_outdir);
    obj.insert("outdir".to_owned(), Value::String(resolved.to_string_lossy().into_owned()));
  }
  let args: S = match serde_json::from_value(body) {
    Ok(args) => args,
    Err(err) => return AppError::from(err).into_response(),
  };
  let real_args: R = args.into();
  sse_response(move |progress| {
    let result = command(&real_args, progress)?;
    serde_json::to_value(result).map_err(Report::from)
  })
}
