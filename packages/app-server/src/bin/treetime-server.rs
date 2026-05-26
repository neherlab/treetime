use app_server::create_router;
use ctor::ctor;
use log::info;
use treetime_utils::init::global::global_init;

#[ctor]
fn init() {
  global_init();
}

#[tokio::main]
async fn main() -> eyre::Result<()> {
  let host = std::env::var("HOST").unwrap_or("127.0.0.1".to_owned());
  let port: u16 = std::env::var("PORT").ok().and_then(|p| p.parse().ok()).unwrap_or(3100);

  let addr = format!("{host}:{port}");
  let listener = tokio::net::TcpListener::bind(&addr).await?;
  info!("TreeTime server listening on http://{addr}");
  axum::serve(listener, create_router()).await?;
  Ok(())
}
