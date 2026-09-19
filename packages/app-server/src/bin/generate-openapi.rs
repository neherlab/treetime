use app_server::routes::api_doc;
use ctor::ctor;
use std::path::PathBuf;
use treetime_utils::init::global::global_init;

#[ctor]
fn init() {
  global_init();
}

fn main() -> eyre::Result<()> {
  let out = std::env::args()
    .nth(1)
    .map_or_else(|| PathBuf::from("packages/app-contracts/openapi.yaml"), PathBuf::from);
  let yaml = api_doc()
    .to_yaml()
    .map_err(|err| eyre::eyre!("failed to serialize the OpenAPI document: {err}"))?;
  std::fs::write(&out, yaml)?;
  eprintln!("Wrote OpenAPI document to {}", out.display());
  Ok(())
}
