use crate::gtr::gtr::GTR;
use crate::nuc_models::jc69::jc69;
use clap::ArgEnum;
use eyre::Report;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
pub enum GtrModelName {
  Infer,
  JC69,
}

impl Default for GtrModelName {
  fn default() -> Self {
    Self::JC69
  }
}

pub fn get_gtr(name: &GtrModelName) -> Result<GTR, Report> {
  match name {
    GtrModelName::Infer => {
      unimplemented!("Not implemented: GTR inference is not yet implemented. Please provide the name of a concrete model with `--gtr=<model>` argument.")
    }
    GtrModelName::JC69 => jc69(),
  }
}
