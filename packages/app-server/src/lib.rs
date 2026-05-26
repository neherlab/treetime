pub mod args;
pub mod progress;
pub mod routes;

use axum::Router;
use tower_http::cors::CorsLayer;
use tower_http::services::{ServeDir, ServeFile};

pub fn create_router() -> Router {
  let api = Router::new().nest("/api", routes::api_routes());

  match std::env::var("STATIC_DIR") {
    Ok(static_dir) => {
      let index = format!("{static_dir}/index.html");
      api.fallback_service(ServeDir::new(&static_dir).fallback(ServeFile::new(index)))
    },
    Err(_) => api,
  }
  .layer(CorsLayer::permissive())
}

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
