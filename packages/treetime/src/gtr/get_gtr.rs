use crate::gtr::gtr::GTR;
use crate::nuc_models::jc69::{f81, hky85, jc69, k80, t92, F81Params, HKY85Params, JC69Params, K80Params, T92Params};
use clap::ArgEnum;
use eyre::Report;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum)]
pub enum GtrModelName {
  Infer,
  F81,
  HKY85,
  JC69,
  K80,
  T92,
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
    GtrModelName::JC69 => jc69(&JC69Params::default()),
    GtrModelName::F81 => f81(&F81Params::default()),
    GtrModelName::HKY85 => hky85(&HKY85Params::default()),
    GtrModelName::K80 => k80(&K80Params::default()),
    GtrModelName::T92 => t92(&T92Params::default()),
  }
}
