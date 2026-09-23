pub mod commands;
pub mod contract;
pub mod error;
pub mod routes;
pub mod sse;
pub mod state;

use crate::state::ServerConfig;
use axum::Router;
use tower_http::cors::CorsLayer;
use tower_http::services::{ServeDir, ServeFile};

pub fn create_router(config: ServerConfig, static_dir: Option<String>) -> Router {
  let api = routes::api_routes(config);

  match static_dir {
    Some(static_dir) => {
      let index = format!("{static_dir}/index.html");
      api.fallback_service(ServeDir::new(&static_dir).fallback(ServeFile::new(index)))
    },
    None => api,
  }
  .layer(CorsLayer::permissive())
}

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor(unsafe)]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
