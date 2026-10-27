use app_backend::bridge::error::ErrorResponse;
use app_server::error::error_response;
use app_server::port::header_pairs;
use axum::Router;
use axum::body::{Body, BodyDataStream, Bytes};
use axum::http::Request;
use axum::response::Response;
use eyre::Report;
use parking_lot::Mutex;
use std::future::Future;
use std::mem;
use std::pin::Pin;
use tokio_stream::StreamExt as _;
use tokio_util::sync::CancellationToken;
use tower::ServiceExt as _;
use treetime_utils::error::report_to_string;
use treetime_utils::{make_error, make_report};
use wasm_bindgen::prelude::{JsError, wasm_bindgen};

#[wasm_bindgen]
pub struct Exchange {
  pub(crate) inner: PortExchange,
}

#[wasm_bindgen]
#[expect(
  tail_expr_drop_order,
  reason = "the exports that wasm-bindgen generates for async methods have the edition of the macro crate; the Rust 2024 migration lint reports their drop order, which this crate does not migrate"
)]
impl Exchange {
  pub async fn head(&self) -> Result<ResponseHead, JsError> {
    self.inner.head().await.map_err(|report| js_error(&report))
  }

  pub async fn chunk(&self) -> Result<Option<Vec<u8>>, JsError> {
    self
      .inner
      .chunk()
      .await
      .map(|chunk| chunk.map(Vec::from))
      .map_err(|report| js_error(&report))
  }

  pub fn abort(&self) {
    self.inner.abort();
  }
}

pub(crate) struct PortExchange {
  state: Mutex<State>,
  cancel: CancellationToken,
}

impl PortExchange {
  pub(crate) fn start(router: Router, request: Result<Request<Body>, Report>) -> Self {
    let response: ResponseFuture = match request {
      Ok(request) => Box::pin(async move {
        match router.oneshot(request).await {
          Ok(response) => response,
          Err(infallible) => match infallible {},
        }
      }),
      Err(report) => Box::pin(async move { error_response(ErrorResponse::from_report(&report)) }),
    };
    Self {
      state: Mutex::new(State::Waiting(response)),
      cancel: CancellationToken::new(),
    }
  }

  pub(crate) async fn head(&self) -> Result<ResponseHead, Report> {
    let response = match self.take()? {
      State::Waiting(response) => response,
      other @ (State::Streaming(_) | State::Reading | State::Ended) => {
        self.replace(other);
        return make_error!("the response head of this exchange was already read");
      },
    };
    let response = self.cancel.run_until_cancelled(response).await.ok_or_else(aborted)?;
    let (head, body) = response.into_parts();
    self.replace(State::Streaming(body.into_data_stream()));
    let (header_names, header_values) = header_pairs(&head.headers).into_iter().unzip();
    Ok(ResponseHead {
      status: head.status.as_u16(),
      header_names,
      header_values,
    })
  }

  pub(crate) async fn chunk(&self) -> Result<Option<Bytes>, Report> {
    let mut body = match self.take()? {
      State::Streaming(body) => body,
      other @ (State::Waiting(_) | State::Reading | State::Ended) => {
        self.replace(other);
        return make_error!("the response head of this exchange must be read before its body");
      },
    };
    match self.cancel.run_until_cancelled(body.next()).await.ok_or_else(aborted)? {
      Some(Ok(bytes)) => {
        self.replace(State::Streaming(body));
        Ok(Some(bytes))
      },
      Some(Err(err)) => {
        self.replace(State::Ended);
        make_error!("When reading the response body: {err}")
      },
      None => {
        self.replace(State::Ended);
        Ok(None)
      },
    }
  }

  pub(crate) fn abort(&self) {
    self.cancel.cancel();
    self.replace(State::Ended);
  }

  fn replace(&self, state: State) -> State {
    mem::replace(&mut *self.state.lock(), state)
  }

  fn take(&self) -> Result<State, Report> {
    if self.cancel.is_cancelled() {
      return Err(aborted());
    }
    match self.replace(State::Reading) {
      State::Reading => make_error!("a read of this exchange is in progress"),
      State::Ended => {
        self.replace(State::Ended);
        make_error!("this exchange has ended")
      },
      state @ (State::Waiting(_) | State::Streaming(_)) => Ok(state),
    }
  }
}

#[wasm_bindgen]
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResponseHead {
  pub status: u16,
  #[wasm_bindgen(getter_with_clone, js_name = headerNames)]
  pub header_names: Vec<String>,
  #[wasm_bindgen(getter_with_clone, js_name = headerValues)]
  pub header_values: Vec<String>,
}

enum State {
  Waiting(ResponseFuture),
  Streaming(BodyDataStream),
  Reading,
  Ended,
}

type ResponseFuture = Pin<Box<dyn Future<Output = Response> + Send>>;

fn aborted() -> Report {
  make_report!("the exchange was aborted")
}

fn js_error(report: &Report) -> JsError {
  JsError::new(&report_to_string(report))
}
