use quick_xml::de::from_reader;
use quick_xml::DeError;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub fn phyloxml_read(reader: impl std::io::Read) -> Result<Phyloxml, DeError> {
  let reader = std::io::BufReader::new(reader);
  from_reader(reader)
}

pub fn phyloxml_write(writer: impl std::io::Write, phyloxml: &Phyloxml) -> std::io::Result<()> {
  details::to_writer_pretty(writer, "phyloxml", phyloxml)
}

mod details {
  use quick_xml::se::to_string_with_root;
  use serde::Serialize;
  use std::io::Cursor;

  /// Adapt `std::fmt::Write` to `std::io::Write`
  pub struct WriteAdapter<W: std::io::Write>(pub W);

  impl<W: std::io::Write> std::fmt::Write for WriteAdapter<W> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
      #[allow(clippy::map_err_ignore)]
      self.0.write_all(s.as_bytes()).map_err(|_| std::fmt::Error)?;
      Ok(())
    }
  }

  pub fn to_writer_pretty<W, T>(writer: W, root_tag: &str, data: &T) -> std::io::Result<()>
  where
    W: std::io::Write,
    T: Serialize,
  {
    let s = to_string_with_root(root_tag, data).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;

    let reader = xml::ParserConfig::new()
      .trim_whitespace(true)
      .ignore_comments(false)
      .create_reader(Cursor::new(s));

    let mut writer = xml::EmitterConfig::new()
      .perform_indent(true)
      .normalize_empty_elements(false)
      .autopad_comments(false)
      .create_writer(writer);

    for event in reader {
      if let Some(event) = event.map_err(to_io)?.as_writer_event() {
        writer.write(event).map_err(to_io)?;
      }
    }
    Ok(())
  }

  fn to_io<E>(e: E) -> std::io::Error
  where
    E: Into<Box<dyn std::error::Error + Send + Sync>>,
  {
    std::io::Error::new(std::io::ErrorKind::Other, e)
  }
}

/// (Incomplete) representation of PhyloXML format, according to schema at http://www.phyloxml.org/1.20/phyloxml.xsd
#[derive(Debug, Serialize, Deserialize)]
pub struct Phyloxml {
  pub phylogeny: Vec<PhyloxmlPhylogeny>,
  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlPhylogeny {
  #[serde(rename = "@rooted")]
  pub rooted: bool,
  #[serde(rename = "@rerootable", skip_serializing_if = "Option::is_none")]
  pub rerootable: Option<bool>,
  #[serde(rename = "@branch_length_unit", skip_serializing_if = "Option::is_none")]
  pub branch_length_unit: Option<String>,
  #[serde(rename = "@type")]
  #[serde(skip_serializing_if = "Option::is_none")]
  pub phylogeny_type: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub name: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub id: Option<PhyloXmlId>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub description: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub date: Option<String>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub confidence: Vec<PhyloxmlConfidence>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub clade: Option<PhyloxmlClade>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub clade_relation: Vec<PhyloxmlCladeRelation>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub sequence_relation: Vec<PhyloxmlSequenceRelation>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub property: Vec<PhyloxmlProperty>,
  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlClade {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub name: Option<String>,
  #[serde(rename = "branch_length")]
  #[serde(skip_serializing_if = "Option::is_none")]
  pub branch_length_elem: Option<f64>,
  #[serde(rename = "@branch_length")]
  #[serde(skip_serializing_if = "Option::is_none")]
  pub branch_length_attr: Option<f64>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub confidence: Vec<PhyloxmlConfidence>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub width: Option<f64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub color: Option<PhyloxmlBranchColor>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub node_id: Option<PhyloXmlId>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub taxonomy: Vec<PhyloxmlTaxonomy>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub sequence: Vec<PhyloxmlSequence>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub events: Option<PhyloxmlEvents>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub binary_characters: Option<PhyloxmlBinaryCharacters>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub distribution: Vec<PhyloxmlDistribution>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub date: Option<PhyloxmlDate>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub reference: Vec<PhyloxmlReference>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub property: Vec<PhyloxmlProperty>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub clade: Vec<PhyloxmlClade>,
  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlTaxonomy {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub id: Option<PhyloXmlId>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub code: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub scientific_name: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub authority: Option<String>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub common_name: Vec<String>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub synonym: Vec<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub rank: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub uri: Option<PhyloxmlUri>,
  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlSequence {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub symbol: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub accession: Option<PhyloxmlAccession>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub name: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub location: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub mol_seq: Option<PhyloxmlMolSeq>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub uri: Option<PhyloxmlUri>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub annotation: Vec<PhyloxmlAnnotation>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub domain_architecture: Option<PhyloxmlDomainArchitecture>,
  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlMolSeq {
  #[serde(rename = "$value")]
  pub sequence: String,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub is_aligned: Option<bool>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlAccession {
  #[serde(rename = "$value")]
  pub accession: String,
  #[serde(rename = "@source")]
  pub source: String,
  #[serde(rename = "@comment")]
  #[serde(skip_serializing_if = "Option::is_none")]
  pub comment: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlDomainArchitecture {
  #[serde(rename = "@length")]
  pub length: u64,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub domain: Vec<PhyloxmlProteinDomain>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlProteinDomain {
  #[serde(rename = "$value")]
  pub name: String,
  #[serde(rename = "@from")]
  pub from: u64,
  #[serde(rename = "@to")]
  pub to: u64,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub confidence: Option<f64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub id: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlEvents {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub event_type: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub duplications: Option<u64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub speciations: Option<u64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub losses: Option<u64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub confidence: Option<PhyloxmlConfidence>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlConfidence {
  #[serde(rename = "$value")]
  pub value: f64,
  #[serde(rename = "@type")]
  pub type_: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloXmlId {
  #[serde(rename = "$value")]
  pub identifier: String,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub provider: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlDistribution {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub desc: Option<String>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub point: Vec<PhyloxmlPoint>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub polygon: Vec<Polygon>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlPoint {
  pub lat: f64,
  pub long: f64,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub alt: Option<f64>,
  pub geodetic_datum: String,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub alt_unit: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Polygon {
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub point: Vec<PhyloxmlPoint>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlDate {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub desc: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub value: Option<f64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub minimum: Option<f64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub maximum: Option<f64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub unit: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlBranchColor {
  pub red: u8,
  pub green: u8,
  pub blue: u8,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlSequenceRelation {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub confidence: Option<PhyloxmlConfidence>,
  pub id_ref_0: String,
  pub id_ref_1: String,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub distance: Option<f64>,
  pub type_: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlCladeRelation {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub confidence: Option<PhyloxmlConfidence>,
  pub id_ref_0: String,
  pub id_ref_1: String,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub distance: Option<f64>,
  pub type_: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlUri {
  #[serde(rename = "$value")]
  pub uri: String,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub desc: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub type_: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlProperty {
  #[serde(rename = "$value")]
  pub value: String,
  #[serde(rename = "@ref")]
  pub ref_: String,
  #[serde(rename = "@unit")]
  #[serde(skip_serializing_if = "Option::is_none")]
  pub unit: Option<String>,
  #[serde(rename = "@datatype")]
  pub datatype: String,
  #[serde(rename = "@applies_to")]
  pub applies_to: String,
  #[serde(rename = "@id_ref")]
  #[serde(skip_serializing_if = "Option::is_none")]
  pub id_ref: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlAnnotation {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub desc: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub confidence: Option<PhyloxmlConfidence>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub property: Vec<PhyloxmlProperty>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub uri: Option<PhyloxmlUri>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub ref_: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub source: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub evidence: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub type_: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlBinaryCharacters {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub gained: Option<PhyloxmlBinaryCharacterList>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub lost: Option<PhyloxmlBinaryCharacterList>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub present: Option<PhyloxmlBinaryCharacterList>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub absent: Option<PhyloxmlBinaryCharacterList>,
  #[serde(rename = "type")]
  #[serde(skip_serializing_if = "Option::is_none")]
  pub character_type: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub gained_count: Option<u64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub lost_count: Option<u64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub present_count: Option<u64>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub absent_count: Option<u64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlBinaryCharacterList {
  #[serde(rename = "bc")]
  pub characters: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlReference {
  #[serde(skip_serializing_if = "Option::is_none")]
  pub desc: Option<String>,
  #[serde(skip_serializing_if = "Option::is_none")]
  pub doi: Option<String>,
}
