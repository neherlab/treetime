use deser::Serialize;
use schemars::JsonSchema;

pub fn version_info() -> VersionInfo {
  VersionInfo {
    version: env!("CARGO_PKG_VERSION"),
  }
}

#[derive(Clone, Debug, JsonSchema, Serialize)]
pub struct VersionInfo {
  pub version: &'static str,
}
