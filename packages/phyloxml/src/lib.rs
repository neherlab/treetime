use quick_xml::de::from_reader;
use quick_xml::se::to_writer;
use quick_xml::DeError;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub fn phyloxml_read(reader: impl std::io::Read) -> Result<Phyloxml, DeError> {
  let reader = std::io::BufReader::new(reader);
  from_reader(reader)
}

pub fn phyloxml_write(writer: impl std::io::Write, phyloxml: &Phyloxml) -> Result<(), DeError> {
  to_writer(&mut details::WriteAdapter(writer), phyloxml)
}

mod details {
  /// Adapt `std::fmt::Write` to `std::io::Write`
  pub struct WriteAdapter<W: std::io::Write>(pub W);

  impl<W: std::io::Write> std::fmt::Write for WriteAdapter<W> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
      #[allow(clippy::map_err_ignore)]
      self.0.write_all(s.as_bytes()).map_err(|_| std::fmt::Error)?;
      Ok(())
    }
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
  #[serde(rename = "@rerootable")]
  pub rerootable: Option<bool>,
  #[serde(rename = "@branch_length_unit")]
  pub branch_length_unit: Option<String>,
  #[serde(rename = "@type")]
  pub phylogeny_type: Option<String>,
  pub name: Option<String>,
  pub id: Option<PhyloXmlId>,
  pub description: Option<String>,
  pub date: Option<String>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub confidence: Vec<PhyloxmlConfidence>,
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
  pub name: Option<String>,
  #[serde(rename = "branch_length")]
  pub branch_length_elem: Option<f64>,
  #[serde(rename = "@branch_length")]
  pub branch_length_attr: Option<f64>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub confidence: Vec<PhyloxmlConfidence>,
  pub width: Option<f64>,
  pub color: Option<PhyloxmlBranchColor>,
  pub node_id: Option<PhyloXmlId>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub taxonomy: Vec<PhyloxmlTaxonomy>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub sequence: Vec<PhyloxmlSequence>,
  pub events: Option<PhyloxmlEvents>,
  pub binary_characters: Option<PhyloxmlBinaryCharacters>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub distribution: Vec<PhyloxmlDistribution>,
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
  pub id: Option<PhyloXmlId>,
  pub code: Option<String>,
  pub scientific_name: Option<String>,
  pub authority: Option<String>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub common_name: Vec<String>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub synonym: Vec<String>,
  pub rank: Option<String>,
  pub uri: Option<PhyloxmlUri>,
  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlSequence {
  pub symbol: Option<String>,
  pub accession: Option<PhyloxmlAccession>,
  pub name: Option<String>,
  pub location: Option<String>,
  pub mol_seq: Option<PhyloxmlMolSeq>,
  pub uri: Option<PhyloxmlUri>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub annotation: Vec<PhyloxmlAnnotation>,
  pub domain_architecture: Option<PhyloxmlDomainArchitecture>,
  #[serde(flatten)]
  pub other: BTreeMap<String, serde_json::Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlMolSeq {
  #[serde(rename = "$value")]
  pub sequence: String,
  pub is_aligned: Option<bool>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlAccession {
  #[serde(rename = "$value")]
  pub accession: String,
  #[serde(rename = "@source")]
  pub source: String,
  #[serde(rename = "@comment")]
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
  pub confidence: Option<f64>,
  pub id: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlEvents {
  pub event_type: Option<String>,
  pub duplications: Option<u64>,
  pub speciations: Option<u64>,
  pub losses: Option<u64>,
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
  pub provider: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlDistribution {
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
  pub alt: Option<f64>,
  pub geodetic_datum: String,
  pub alt_unit: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Polygon {
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub point: Vec<PhyloxmlPoint>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlDate {
  pub desc: Option<String>,
  pub value: Option<f64>,
  pub minimum: Option<f64>,
  pub maximum: Option<f64>,
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
  pub confidence: Option<PhyloxmlConfidence>,
  pub id_ref_0: String,
  pub id_ref_1: String,
  pub distance: Option<f64>,
  pub type_: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlCladeRelation {
  pub confidence: Option<PhyloxmlConfidence>,
  pub id_ref_0: String,
  pub id_ref_1: String,
  pub distance: Option<f64>,
  pub type_: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlUri {
  #[serde(rename = "$value")]
  pub uri: String,
  pub desc: Option<String>,
  pub type_: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlProperty {
  #[serde(rename = "$value")]
  pub value: String,
  #[serde(rename = "@ref")]
  pub ref_: String,
  #[serde(rename = "@unit")]
  pub unit: Option<String>,
  #[serde(rename = "@datatype")]
  pub datatype: String,
  #[serde(rename = "@applies_to")]
  pub applies_to: String,
  #[serde(rename = "@id_ref")]
  pub id_ref: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlAnnotation {
  pub desc: Option<String>,
  pub confidence: Option<PhyloxmlConfidence>,
  #[serde(default, skip_serializing_if = "Vec::is_empty")]
  pub property: Vec<PhyloxmlProperty>,
  pub uri: Option<PhyloxmlUri>,
  pub ref_: Option<String>,
  pub source: Option<String>,
  pub evidence: Option<String>,
  pub type_: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlBinaryCharacters {
  pub gained: Option<PhyloxmlBinaryCharacterList>,
  pub lost: Option<PhyloxmlBinaryCharacterList>,
  pub present: Option<PhyloxmlBinaryCharacterList>,
  pub absent: Option<PhyloxmlBinaryCharacterList>,
  #[serde(rename = "type")]
  pub character_type: Option<String>,
  pub gained_count: Option<u64>,
  pub lost_count: Option<u64>,
  pub present_count: Option<u64>,
  pub absent_count: Option<u64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlBinaryCharacterList {
  #[serde(rename = "bc")]
  pub characters: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlReference {
  pub desc: Option<String>,
  pub doi: Option<String>,
}
