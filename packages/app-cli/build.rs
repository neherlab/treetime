use std::env;
use std::error::Error;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use treetime_schema::{TreetimeSchemaFormat, generate_schema};

fn main() -> Result<(), Box<dyn Error>> {
  println!("cargo:rerun-if-changed=../treetime/src");

  let out_dir = PathBuf::from("../app-contracts/src/generated");
  if let Err(err) = generate_schema(&TreetimeSchemaFormat::All, Some(&out_dir)) {
    eprintln!("cargo:warning=Schema generation failed: {err}");
  }

  emit_long_version()?;
  Ok(())
}

fn emit_long_version() -> Result<(), Box<dyn Error>> {
  println!("cargo:rerun-if-env-changed=TREETIME_BUILD_MODE");
  println!("cargo:rerun-if-env-changed=TREETIME_VERSION_SUFFIX");
  emit_git_dependency("HEAD")?;
  emit_git_dependency("index")?;

  let symbolic_ref = git_output(&["symbolic-ref", "--quiet", "HEAD"])?;
  if symbolic_ref.status.success() {
    emit_git_dependency(&output_stdout(&symbolic_ref)?)?;
  } else if symbolic_ref.status.code() != Some(1) {
    return Err(command_error("git symbolic-ref --quiet HEAD", &symbolic_ref).into());
  }

  let base = env!("CARGO_PKG_VERSION");
  let mode = env::var("TREETIME_BUILD_MODE").unwrap_or_else(|_| "dev".to_owned());
  let suffix = env::var("TREETIME_VERSION_SUFFIX")
    .ok()
    .filter(|suffix| !suffix.is_empty());
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

fn git_output(args: &[&str]) -> Result<Output, std::io::Error> {
  Command::new("git").args(args).output()
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
