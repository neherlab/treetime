use crate::commands::homoplasy::result::DrmAnnotation;
use eyre::{Report, WrapErr};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::Path;
use treetime::make_error;
use treetime_io::csv::{TableFormat, csv_read_file};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct DrmTable {
  positions: BTreeMap<usize, DrmPosition>,
}

impl DrmTable {
  pub fn read_file(path: &Path) -> Result<Self, Report> {
    csv_read_file(path, TableFormat::Tsv)
      .and_then(Self::from_rows)
      .wrap_err_with(|| format!("When reading drug resistance mutations from '{}'", path.display()))
  }

  pub fn from_rows(rows: Vec<DrmRow>) -> Result<Self, Report> {
    let mut positions: BTreeMap<usize, DrmPosition> = BTreeMap::new();
    for row in rows {
      let Some(position) = row.genomic_position.checked_sub(1) else {
        return make_error!("GENOMIC_POSITION counts from 1, but a row has position 0");
      };
      positions
        .entry(position)
        .or_insert_with(|| DrmPosition {
          gene: row.gene,
          drug: row.drug,
          substitutions: BTreeMap::new(),
        })
        .substitutions
        .entry(row.alt_base)
        .or_insert(row.substitution);
    }
    Ok(Self { positions })
  }

  pub fn contains(&self, position: usize) -> bool {
    self.positions.contains_key(&position)
  }

  pub fn annotate(&self, position: usize, derived: &str) -> Option<DrmAnnotation> {
    self.positions.get(&position).map(|drm| DrmAnnotation {
      gene: drm.gene.clone(),
      drug: drm.drug.clone(),
      substitution: drm.substitutions.get(derived).cloned(),
    })
  }
}

#[derive(Clone, Debug, PartialEq, Eq, Deserialize)]
pub struct DrmRow {
  #[serde(rename = "GENOMIC_POSITION")]
  pub genomic_position: usize,
  #[serde(rename = "ALT_BASE")]
  pub alt_base: String,
  #[serde(rename = "DRUG")]
  pub drug: String,
  #[serde(rename = "GENE")]
  pub gene: String,
  #[serde(rename = "SUBSTITUTION")]
  pub substitution: String,
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct DrmPosition {
  gene: String,
  drug: String,
  substitutions: BTreeMap<String, String>,
}
