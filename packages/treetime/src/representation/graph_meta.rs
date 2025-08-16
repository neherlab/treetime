use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MetaNodePartition;

impl MetaNodePartition {
  pub fn new() -> Result<Self, eyre::Report> {
    Ok(Self {})
  }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct MetaEdgePartition;

impl MetaEdgePartition {
  pub fn new() -> Result<Self, eyre::Report> {
    Ok(Self {})
  }
}
