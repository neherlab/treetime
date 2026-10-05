use app_commands::config::catalog::setting_catalog;
use ctor::ctor;
use std::path::PathBuf;
use treetime_utils::init::global::global_init;
use treetime_utils::io::json::{JsonPretty, json_write_file};

#[ctor(unsafe)]
fn init() {
  global_init();
}

fn main() -> eyre::Result<()> {
  let out = std::env::args().nth(1).map_or_else(
    || PathBuf::from("packages/app-contracts/setting-catalog.json"),
    PathBuf::from,
  );
  json_write_file(&out, &setting_catalog()?, JsonPretty(true))
}
