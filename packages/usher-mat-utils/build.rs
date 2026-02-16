use eyre::Report;

fn main() -> Result<(), Report> {
  prost_build::Config::new()
    .type_attribute(".", "#[must_use]")
    .type_attribute(".", "#[derive(serde::Serialize,serde::Deserialize)]")
    .compile_protos(&["schemas/parsimony.proto"], &["src/gen/"])?;
  Ok(())
}
