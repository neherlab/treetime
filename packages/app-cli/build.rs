#![cfg_attr(
  dylint_lib = "treetime_lints",
  expect(debug_remnants, reason = "cargo reads build script directives from stdout")
)]

use std::env::{self, VarError};
use std::error::Error;
use std::path::Path;
use std::process::{Command, Output};

const BUILD_MODE_ENV: &str = "TREETIME_BUILD_MODE";
const VERSION_SUFFIX_ENV: &str = "TREETIME_VERSION_SUFFIX";

fn main() -> Result<(), Box<dyn Error>> {
  emit_long_version()?;
  Ok(())
}

fn emit_long_version() -> Result<(), Box<dyn Error>> {
  println!("cargo:rerun-if-env-changed={BUILD_MODE_ENV}");
  println!("cargo:rerun-if-env-changed={VERSION_SUFFIX_ENV}");
  emit_git_dependency("HEAD")?;
  emit_git_dependency("index")?;

  let symbolic_ref = git_output(&["symbolic-ref", "--quiet", "HEAD"])?;
  if symbolic_ref.status.success() {
    emit_git_dependency(&output_stdout(&symbolic_ref)?)?;
  } else if symbolic_ref.status.code() != Some(1) {
    return Err(command_error("git symbolic-ref --quiet HEAD", &symbolic_ref).into());
  }

  let base = env!("CARGO_PKG_VERSION");
  let mode = env_var_optional(BUILD_MODE_ENV)?.unwrap_or_else(|| "dev".to_owned());
  let suffix = env_var_optional(VERSION_SUFFIX_ENV)?.filter(|suffix| !suffix.is_empty());
  let short_sha = git_stdout(&["rev-parse", "--short", "HEAD"])?;
  let dirty = if git_stdout(&["status", "--porcelain"])?.is_empty() {
    ""
  } else {
    ".dirty"
  };

  let long_version = match (mode.as_str(), suffix) {
    ("dev", None) => format!("{base}-dev+{short_sha}{dirty}"),
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
  Ok(())
}

fn env_var_optional(name: &str) -> Result<Option<String>, VarError> {
  match env::var(name) {
    Ok(value) => Ok(Some(value)),
    Err(VarError::NotPresent) => Ok(None),
    Err(error @ VarError::NotUnicode(_)) => Err(error),
  }
}

fn emit_git_dependency(path: &str) -> Result<(), Box<dyn Error>> {
  let path = git_stdout(&["rev-parse", "--path-format=absolute", "--git-path", path])?;
  if Path::new(&path).exists() {
    println!("cargo:rerun-if-changed={path}");
  }
  Ok(())
}

fn git_stdout(args: &[&str]) -> Result<String, Box<dyn Error>> {
  let output = git_output(args)?;
  if !output.status.success() {
    return Err(command_error(&format!("git {}", args.join(" ")), &output).into());
  }
  Ok(output_stdout(&output)?)
}

fn git_output(args: &[&str]) -> Result<Output, Box<dyn Error>> {
  Command::new("git")
    .args(args)
    .output()
    .map_err(|error| format!("When running git {}: {error}", args.join(" ")).into())
}

fn output_stdout(output: &Output) -> Result<String, std::str::Utf8Error> {
  Ok(std::str::from_utf8(&output.stdout)?.trim().to_owned())
}

fn command_error(command: &str, output: &Output) -> String {
  format!(
    "{command} failed with {}: {}",
    output.status,
    String::from_utf8_lossy(&output.stderr).trim()
  )
}
