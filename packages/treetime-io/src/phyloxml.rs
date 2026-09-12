use eyre::{Report, WrapErr};
use smart_default::SmartDefault;
use std::io::Write;
use std::path::Path;
use treetime_utils::io::file::create_file_or_stdout;
use treetime_utils::io::json::{JsonPretty, json_write, json_write_file, json_write_str};
pub use util_phyloxml::{
  Phyloxml, PhyloxmlAccession, PhyloxmlAnnotation, PhyloxmlBinaryCharacterList, PhyloxmlBinaryCharacters,
  PhyloxmlBranchColor, PhyloxmlClade, PhyloxmlCladeRelation, PhyloxmlConfidence, PhyloxmlDate, PhyloxmlDistribution,
  PhyloxmlDomainArchitecture, PhyloxmlEvents, PhyloxmlId, PhyloxmlMolSeq, PhyloxmlPhylogeny, PhyloxmlPoint,
  PhyloxmlProperty, PhyloxmlProteinDomain, PhyloxmlReference, PhyloxmlSequence, PhyloxmlSequenceRelation,
  PhyloxmlTaxonomy, PhyloxmlUri, Polygon,
};

pub fn phyloxml_write_file(filepath: impl AsRef<Path>, phyloxml: &Phyloxml) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  let mut f = create_file_or_stdout(filepath)?;
  phyloxml_write(&mut f, phyloxml).wrap_err_with(|| format!("When writing PhyloXML file '{}'", filepath.display()))?;
  writeln!(f)?;
  Ok(())
}

pub fn phyloxml_write_str(phyloxml: &Phyloxml) -> Result<String, Report> {
  let mut buf = Vec::new();
  phyloxml_write(&mut buf, phyloxml).wrap_err("When writing PhyloXML string")?;
  String::from_utf8(buf).wrap_err("PhyloXML output is not valid UTF-8")
}

pub fn phyloxml_write(writer: &mut impl Write, phyloxml: &Phyloxml) -> Result<(), Report> {
  util_phyloxml::phyloxml_write(writer, phyloxml).wrap_err("When writing PhyloXML")
}

pub fn phyloxml_json_write_file(
  filepath: impl AsRef<Path>,
  phyloxml: &Phyloxml,
  options: &PhyloxmlJsonOptions,
) -> Result<(), Report> {
  let filepath = filepath.as_ref();
  json_write_file(filepath, phyloxml, JsonPretty(options.pretty))
    .wrap_err_with(|| format!("When writing PhyloXML JSON file: '{}'", filepath.display()))?;
  Ok(())
}

pub fn phyloxml_json_write_str(phyloxml: &Phyloxml, options: &PhyloxmlJsonOptions) -> Result<String, Report> {
  json_write_str(phyloxml, JsonPretty(options.pretty)).wrap_err("When writing PhyloXML JSON string")
}

pub fn phyloxml_json_write(
  writer: &mut impl Write,
  phyloxml: &Phyloxml,
  options: &PhyloxmlJsonOptions,
) -> Result<(), Report> {
  json_write(writer, phyloxml, JsonPretty(options.pretty)).wrap_err("When writing PhyloXML JSON")
}

#[derive(SmartDefault)]
pub struct PhyloxmlJsonOptions {
  #[default = true]
  pretty: bool,
}

pub struct PhyloxmlNodeImpl {
  pub name: Option<String>,
  pub branch_length: f64,
}
