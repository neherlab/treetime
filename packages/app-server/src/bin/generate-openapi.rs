use app_server::routes::api_doc;
use ctor::ctor;
use eyre::WrapErr;
use std::fs;
use std::path::PathBuf;
use treetime_utils::init::global::global_init;
use treetime_utils::make_report;

#[ctor(unsafe)]
fn init() {
  global_init();
}

fn main() -> eyre::Result<()> {
  let out = std::env::args()
    .nth(1)
    .map_or_else(|| PathBuf::from("packages/app-contracts/openapi.yaml"), PathBuf::from);
  let yaml = api_doc()
    .to_yaml()
    .map_err(|err| make_report!("When serializing the OpenAPI document: {err}"))?;
  fs::write(&out, yaml).wrap_err_with(|| format!("When writing the OpenAPI document to '{}'", out.display()))?;
  Ok(())
}
