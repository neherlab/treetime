#![cfg_attr(
  dylint_lib = "treetime_lints",
  expect(debug_remnants, reason = "cargo reads build script directives from stdout")
)]

use std::env::{self, VarError};
use std::error::Error;

const BUILD_MODE_ENV: &str = "TREETIME_BUILD_MODE";
const VERSION_SUFFIX_ENV: &str = "TREETIME_VERSION_SUFFIX";

fn main() -> Result<(), Box<dyn Error>> {
  println!("cargo:rerun-if-changed=build.rs");
  println!("cargo:rerun-if-env-changed={BUILD_MODE_ENV}");
  println!("cargo:rerun-if-env-changed={VERSION_SUFFIX_ENV}");

  let base = env!("CARGO_PKG_VERSION");
  let mode = env_var_optional(BUILD_MODE_ENV)?.unwrap_or_else(|| "dev".to_owned());
  let suffix = env_var_optional(VERSION_SUFFIX_ENV)?.filter(|suffix| !suffix.is_empty());

  let long_version = match (mode.as_str(), suffix) {
    ("dev", None) => format!("{base}-dev"),
    ("nightly", Some(suffix)) => format!("{base}-{suffix}"),
    ("release", None) => base.to_owned(),
    ("dev" | "release", Some(_)) => {
      return Err(format!("TREETIME_VERSION_SUFFIX is invalid for {mode} builds").into());
    },
    ("nightly", None) => {
      return Err("TREETIME_VERSION_SUFFIX is required for nightly builds".into());
    },
    _ => return Err(format!("Unknown TREETIME_BUILD_MODE: {mode}").into()),
  };

  println!("cargo:rustc-env=TREETIME_LONG_VERSION={long_version}");
  println!("cargo:rustc-env=TREETIME_BUILD_MODE={mode}");
  Ok(())
}

fn env_var_optional(name: &str) -> Result<Option<String>, VarError> {
  match env::var(name) {
    Ok(value) => Ok(Some(value)),
    Err(VarError::NotPresent) => Ok(None),
    Err(error @ VarError::NotUnicode(_)) => Err(error),
  }
}
