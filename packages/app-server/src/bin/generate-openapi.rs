use app_server::routes::api_doc;
use ctor::ctor;
use deser::adapters::As;
use deser_serde::Serde;
use std::path::PathBuf;
use treetime_utils::init::global::global_init;
use treetime_utils::io::json::{JsonPretty, json_write_file};

#[ctor(unsafe)]
fn init() {
  global_init();
}

fn main() -> eyre::Result<()> {
  let out = std::env::args()
    .nth(1)
    .map_or_else(|| PathBuf::from("packages/app-contracts/openapi.json"), PathBuf::from);
  json_write_file(&out, &As::<_, Serde>::new(&api_doc()?), JsonPretty(true))
}
