#[allow(clippy::needless_pass_by_value)]
pub mod routes;

use axum::Router;
use tower_http::cors::CorsLayer;

pub fn create_router() -> Router {
  Router::new()
    .nest("/api", routes::api_routes())
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
