use eyre::{Report, WrapErr};

fn main() -> Result<(), Report> {
  prost_build::Config::new()
    .type_attribute(".", "#[must_use]")
    .type_attribute(".", "#[derive(serde::Serialize,serde::Deserialize)]")
    .compile_protos(
      &[
        "schemas/mutation_detailed.proto",
        "schemas/parsimony.proto",
        "schemas/taxodium.proto",
      ],
      &["schemas/"],
    )
    .wrap_err("When compiling the UShER protobuf schemas")?;
  Ok(())
}
