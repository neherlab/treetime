use crate::alphabet::alphabet::Alphabet;
use crate::representation::graph_dense::DenseNodePartition;
use crate::representation::graph_meta::MetaNodePartition;
use crate::representation::graph_sparse::SparseNodePartition;
use crate::representation::seq::Seq;
use enum_extract_macro::EnumExtract;
use eyre::Report;
use serde::{Deserialize, Serialize};
use strum_macros::IntoStaticStr;

#[derive(Clone, Debug, Serialize, Deserialize, IntoStaticStr, EnumExtract)]
pub enum NodePartition {
  Dense(DenseNodePartition),
  Sparse(SparseNodePartition),
  Meta(MetaNodePartition),
}

impl NodePartition {
  pub fn dense(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
    Ok(NodePartition::Dense(DenseNodePartition::new(seq, alphabet)?))
  }

  pub fn sparse(seq: &Seq, alphabet: &Alphabet) -> Result<Self, Report> {
    Ok(NodePartition::Sparse(SparseNodePartition::new(seq, alphabet)?))
  }

  pub fn meta() -> Result<Self, Report> {
    Ok(NodePartition::Meta(MetaNodePartition::new()?))
  }

  // pub fn is_dense(&self) -> bool {
  //   matches!(self, NodePartition::Dense(_))
  // }

  // pub fn is_sparse(&self) -> bool {
  //   matches!(self, NodePartition::Sparse(_))
  // }

  // pub fn is_meta(&self) -> bool {
  //   matches!(self, NodePartition::Meta(_))
  // }

  // pub fn as_dense(&self) -> Result<&DenseNodePartition, Report> {
  //   match self {
  //     NodePartition::Dense(partition) => Ok(partition),
  //     _ => make_internal_error!("Expected NodePartition to be Dense, but got {}", self.variant_name()),
  //   }
  // }

  // pub fn as_sparse(&self) -> Result<&SparseNodePartition, Report> {
  //   match self {
  //     NodePartition::Sparse(partition) => Ok(partition),
  //     _ => make_internal_error!("Expected NodePartition to be Sparse, but got {}", self.variant_name()),
  //   }
  // }

  // pub fn as_meta(&self) -> Result<&MetaNodePartition, Report> {
  //   match self {
  //     NodePartition::Meta(partition) => Ok(partition),
  //     _ => make_internal_error!("Expected NodePartition to be Meta, but got {}", self.variant_name()),
  //   }
  // }

  // pub fn as_dense_mut(&mut self) -> Result<&mut DenseNodePartition, Report> {
  //   match self {
  //     NodePartition::Dense(partition) => Ok(partition),
  //     _ => make_internal_error!("Expected NodePartition to be Dense, but got {}", self.variant_name()),
  //   }
  // }

  // pub fn as_sparse_mut(&mut self) -> Result<&mut SparseNodePartition, Report> {
  //   match self {
  //     NodePartition::Sparse(partition) => Ok(partition),
  //     _ => make_internal_error!("Expected NodePartition to be Sparse, but got {}", self.variant_name()),
  //   }
  // }

  // pub fn as_meta_mut(&mut self) -> Result<&mut MetaNodePartition, Report> {
  //   match self {
  //     NodePartition::Meta(partition) => Ok(partition),
  //     _ => make_internal_error!("Expected NodePartition to be Meta, but got {}", self.variant_name()),
  //   }
  // }

  fn variant_name(&self) -> &'static str {
    self.into()
  }
}
