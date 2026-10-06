use schemars::JsonSchema;
use serde::Serialize;

pub fn version_info() -> VersionInfo {
  VersionInfo {
    version: env!("CARGO_PKG_VERSION"),
  }
}

#[derive(Clone, Debug, Serialize, JsonSchema, deser::Serialize)]
pub struct VersionInfo {
  pub version: &'static str,
}
