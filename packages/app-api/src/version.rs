use serde::Serialize;

#[derive(Clone, Debug, Serialize)]
pub struct VersionInfo {
  pub version: &'static str,
}

pub fn version_info() -> VersionInfo {
  VersionInfo {
    version: env!("CARGO_PKG_VERSION"),
  }
}
