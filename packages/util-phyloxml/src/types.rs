use deser::{Deserialize, Serialize};
use deser_value::Value;
use std::collections::BTreeMap;

#[derive(Debug, Serialize, Deserialize)]
#[deser(rename = "phyloxml")]
pub struct Phyloxml {
  pub phylogeny: Vec<PhyloxmlPhylogeny>,
  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlPhylogeny {
  #[deser(rename = "@rooted")]
  pub rooted: bool,
  #[deser(rename = "@rerootable", skip_serializing_if = Option::is_none)]
  pub rerootable: Option<bool>,
  #[deser(rename = "@branch_length_unit", skip_serializing_if = Option::is_none)]
  pub branch_length_unit: Option<String>,
  #[deser(rename = "@type")]
  #[deser(skip_serializing_if = Option::is_none)]
  pub phylogeny_type: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub name: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub id: Option<PhyloxmlId>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub description: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub date: Option<String>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub confidence: Vec<PhyloxmlConfidence>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub clade: Option<PhyloxmlClade>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub clade_relation: Vec<PhyloxmlCladeRelation>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub sequence_relation: Vec<PhyloxmlSequenceRelation>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub property: Vec<PhyloxmlProperty>,
  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlClade {
  #[deser(skip_serializing_if = Option::is_none)]
  pub name: Option<String>,
  #[deser(rename = "branch_length")]
  #[deser(skip_serializing_if = Option::is_none)]
  pub branch_length_elem: Option<f64>,
  #[deser(rename = "@branch_length")]
  #[deser(skip_serializing_if = Option::is_none)]
  pub branch_length_attr: Option<f64>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub confidence: Vec<PhyloxmlConfidence>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub width: Option<f64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub color: Option<PhyloxmlBranchColor>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub node_id: Option<PhyloxmlId>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub taxonomy: Vec<PhyloxmlTaxonomy>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub sequence: Vec<PhyloxmlSequence>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub events: Option<PhyloxmlEvents>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub binary_characters: Option<PhyloxmlBinaryCharacters>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub distribution: Vec<PhyloxmlDistribution>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub date: Option<PhyloxmlDate>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub reference: Vec<PhyloxmlReference>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub property: Vec<PhyloxmlProperty>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub clade: Vec<PhyloxmlClade>,
  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlTaxonomy {
  #[deser(skip_serializing_if = Option::is_none)]
  pub id: Option<PhyloxmlId>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub code: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub scientific_name: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub authority: Option<String>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub common_name: Vec<String>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub synonym: Vec<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub rank: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub uri: Option<PhyloxmlUri>,
  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlSequence {
  #[deser(skip_serializing_if = Option::is_none)]
  pub symbol: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub accession: Option<PhyloxmlAccession>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub name: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub location: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub mol_seq: Option<PhyloxmlMolSeq>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub uri: Option<PhyloxmlUri>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub annotation: Vec<PhyloxmlAnnotation>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub domain_architecture: Option<PhyloxmlDomainArchitecture>,
  #[deser(flatten)]
  pub other: BTreeMap<String, Value>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlMolSeq {
  #[deser(rename = "$text")]
  pub sequence: String,
  #[deser(skip_serializing_if = Option::is_none)]
  pub is_aligned: Option<bool>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlAccession {
  #[deser(rename = "$text")]
  pub accession: String,
  #[deser(rename = "@source")]
  pub source: String,
  #[deser(rename = "@comment")]
  #[deser(skip_serializing_if = Option::is_none)]
  pub comment: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlDomainArchitecture {
  #[deser(rename = "@length")]
  pub length: u64,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub domain: Vec<PhyloxmlProteinDomain>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlProteinDomain {
  #[deser(rename = "$text")]
  pub name: String,
  #[deser(rename = "@from")]
  pub from: u64,
  #[deser(rename = "@to")]
  pub to: u64,
  #[deser(skip_serializing_if = Option::is_none)]
  pub confidence: Option<f64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub id: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlEvents {
  #[deser(skip_serializing_if = Option::is_none)]
  pub event_type: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub duplications: Option<u64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub speciations: Option<u64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub losses: Option<u64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub confidence: Option<PhyloxmlConfidence>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlId {
  #[deser(rename = "$text")]
  pub identifier: String,
  #[deser(skip_serializing_if = Option::is_none)]
  pub provider: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlDistribution {
  #[deser(skip_serializing_if = Option::is_none)]
  pub desc: Option<String>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub point: Vec<PhyloxmlPoint>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub polygon: Vec<Polygon>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct Polygon {
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub point: Vec<PhyloxmlPoint>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlPoint {
  pub lat: f64,
  pub long: f64,
  #[deser(skip_serializing_if = Option::is_none)]
  pub alt: Option<f64>,
  pub geodetic_datum: String,
  #[deser(skip_serializing_if = Option::is_none)]
  pub alt_unit: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlDate {
  #[deser(skip_serializing_if = Option::is_none)]
  pub desc: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub value: Option<f64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub minimum: Option<f64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub maximum: Option<f64>,
  #[deser(skip_serializing_if = Option::is_none)]
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
  #[deser(skip_serializing_if = Option::is_none)]
  pub confidence: Option<PhyloxmlConfidence>,
  pub id_ref_0: String,
  pub id_ref_1: String,
  #[deser(skip_serializing_if = Option::is_none)]
  pub distance: Option<f64>,
  pub type_: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlCladeRelation {
  #[deser(skip_serializing_if = Option::is_none)]
  pub confidence: Option<PhyloxmlConfidence>,
  pub id_ref_0: String,
  pub id_ref_1: String,
  #[deser(skip_serializing_if = Option::is_none)]
  pub distance: Option<f64>,
  pub type_: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlAnnotation {
  #[deser(skip_serializing_if = Option::is_none)]
  pub desc: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub confidence: Option<PhyloxmlConfidence>,
  #[deser(default, skip_serializing_if = Vec::is_empty)]
  pub property: Vec<PhyloxmlProperty>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub uri: Option<PhyloxmlUri>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub ref_: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub source: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub evidence: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub type_: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlConfidence {
  #[deser(rename = "$text")]
  pub value: f64,
  #[deser(rename = "@type")]
  pub type_: String,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlUri {
  #[deser(rename = "$text")]
  pub uri: String,
  #[deser(skip_serializing_if = Option::is_none)]
  pub desc: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub type_: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlProperty {
  #[deser(rename = "$text")]
  pub value: String,
  #[deser(rename = "@ref")]
  pub ref_: String,
  #[deser(rename = "@unit")]
  #[deser(skip_serializing_if = Option::is_none)]
  pub unit: Option<String>,
  #[deser(rename = "@datatype")]
  pub datatype: String,
  #[deser(rename = "@applies_to")]
  pub applies_to: String,
  #[deser(rename = "@id_ref")]
  #[deser(skip_serializing_if = Option::is_none)]
  pub id_ref: Option<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlBinaryCharacters {
  #[deser(skip_serializing_if = Option::is_none)]
  pub gained: Option<PhyloxmlBinaryCharacterList>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub lost: Option<PhyloxmlBinaryCharacterList>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub present: Option<PhyloxmlBinaryCharacterList>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub absent: Option<PhyloxmlBinaryCharacterList>,
  #[deser(rename = "type")]
  #[deser(skip_serializing_if = Option::is_none)]
  pub character_type: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub gained_count: Option<u64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub lost_count: Option<u64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub present_count: Option<u64>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub absent_count: Option<u64>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlBinaryCharacterList {
  #[deser(rename = "bc")]
  pub characters: Vec<String>,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct PhyloxmlReference {
  #[deser(skip_serializing_if = Option::is_none)]
  pub desc: Option<String>,
  #[deser(skip_serializing_if = Option::is_none)]
  pub doi: Option<String>,
}
