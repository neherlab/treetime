use std::path::PathBuf;
use treetime::schema::{TreetimeSchemaFormat, generate_schema};

fn main() {
  println!("cargo:rerun-if-changed=../treetime/src");

  let out_dir = PathBuf::from("../app-contracts/src/generated");
  if let Err(err) = generate_schema(&TreetimeSchemaFormat::All, Some(&out_dir)) {
    eprintln!("cargo:warning=Schema generation failed: {err}");
  }
}
