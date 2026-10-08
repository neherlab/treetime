#[cfg(test)]
mod tests {
  use crate::exchange::{PortExchange, ResponseHead};
  use app_commands::runs::errors::invalid;
  use axum::Router;
  use axum::body::{Body, Bytes};
  use axum::http::Request;
  use axum::routing::get;
  use helpers::{exchange, read_all, test_router};
  use pretty_assertions::assert_eq;
  use serde_json::{Value, json};
  use std::io;
  use tokio_stream::{iter, pending};
  use treetime_utils::{assert_error, o, vec_of_owned};

  #[tokio::test]
  async fn test_exchange_delivers_the_head_then_each_chunk_then_the_end() {
    let exchange = exchange(test_router(), "/chunks");
    let expected = (
      ResponseHead {
        status: 200,
        header_names: vec_of_owned!["cache-control"],
        header_values: vec_of_owned!["no-store"],
      },
      vec![Bytes::from("one"), Bytes::from("two")],
    );
    assert_eq!(expected, read_all(&exchange).await);
  }

  #[tokio::test]
  async fn test_exchange_answers_a_request_that_cannot_be_built_with_an_invalid_request_response() {
    let exchange = PortExchange::start(test_router(), Err(invalid("the request is malformed")));
    let (head, chunks) = read_all(&exchange).await;
    let body: Value = serde_json::from_slice(&chunks.concat()).unwrap();
    let expected = (
      400,
      json!({ "code": "invalid_request", "message": "the request is malformed", "causes": [] }),
    );
    assert_eq!(expected, (head.status, body));
  }

  #[tokio::test]
  async fn test_exchange_refuses_the_body_before_the_head() {
    let exchange = exchange(test_router(), "/chunks");
    assert_error!(
      exchange.chunk().await,
      "the response head of this exchange must be read before its body"
    );
  }

  #[tokio::test]
  async fn test_exchange_refuses_a_second_head() {
    let exchange = exchange(test_router(), "/chunks");
    exchange.head().await.unwrap();
    assert_error!(
      exchange.head().await,
      "the response head of this exchange was already read"
    );
  }

  #[tokio::test]
  async fn test_exchange_refuses_reads_after_the_end() {
    let exchange = exchange(test_router(), "/chunks");
    read_all(&exchange).await;
    assert_error!(exchange.chunk().await, "this exchange has ended");
  }

  #[tokio::test]
  async fn test_exchange_refuses_reads_after_an_abort() {
    let exchange = exchange(test_router(), "/chunks");
    exchange.abort();
    assert_error!(exchange.head().await, "the exchange was aborted");
  }

  #[tokio::test]
  async fn test_exchange_abort_ends_a_read_that_waits_for_a_chunk() {
    let exchange = exchange(test_router(), "/endless");
    exchange.head().await.unwrap();
    let (read, ()) = tokio::join!(exchange.chunk(), async { exchange.abort() });
    assert_error!(read, "the exchange was aborted");
  }

  #[tokio::test]
  async fn test_exchange_refuses_a_second_read_while_one_waits() {
    let exchange = exchange(test_router(), "/endless");
    exchange.head().await.unwrap();
    let (_, second, ()) = tokio::join!(exchange.chunk(), exchange.chunk(), async { exchange.abort() });
    assert_error!(second, "a read of this exchange is in progress");
  }

  #[tokio::test]
  async fn test_exchange_body_failure_ends_the_exchange_with_the_error() {
    let exchange = exchange(test_router(), "/broken");
    exchange.head().await.unwrap();
    let failed = exchange.chunk().await;
    let after = exchange.chunk().await;
    assert_error!(failed, "When reading the response body: the body broke");
    assert_error!(after, "this exchange has ended");
  }

  mod helpers {
    use super::{Body, Bytes, PortExchange, Request, ResponseHead, Router, get, io, iter, o, pending};
    use axum::response::IntoResponse;

    pub(super) fn test_router() -> Router {
      Router::new()
        .route(
          "/chunks",
          get(|| async {
            (
              [("cache-control", "no-store")],
              Body::from_stream(iter([Ok::<_, io::Error>(o!("one")), Ok(o!("two"))])),
            )
              .into_response()
          }),
        )
        .route(
          "/endless",
          get(|| async { Body::from_stream(pending::<Result<Bytes, io::Error>>()) }),
        )
        .route(
          "/broken",
          get(|| async { Body::from_stream(iter([Err::<Bytes, _>(io::Error::other("the body broke"))])) }),
        )
    }

    pub(super) fn exchange(router: Router, url: &str) -> PortExchange {
      PortExchange::start(router, Ok(Request::get(url).body(Body::empty()).unwrap()))
    }

    pub(super) async fn read_all(exchange: &PortExchange) -> (ResponseHead, Vec<Bytes>) {
      let head = exchange.head().await.unwrap();
      let mut chunks = vec![];
      while let Some(chunk) = exchange.chunk().await.unwrap() {
        chunks.push(chunk);
      }
      (head, chunks)
    }
  }
}
