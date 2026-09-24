use app_server::routes::api_doc;
use ctor::ctor;
use eyre::WrapErr;
use std::fs;
use std::path::PathBuf;
use treetime_utils::init::global::global_init;
use treetime_utils::io::json::{JsonPretty, json_write_str};

#[ctor(unsafe)]
fn init() {
  global_init();
}

fn main() -> eyre::Result<()> {
  let out = std::env::args()
    .nth(1)
    .map_or_else(|| PathBuf::from("packages/app-contracts/openapi.json"), PathBuf::from);
  let doc = json_write_str(&api_doc()?, JsonPretty(true))?;
  fs::write(&out, format!("{doc}\n"))
    .wrap_err_with(|| format!("When writing the OpenAPI document to '{}'", out.display()))?;
  Ok(())
}
