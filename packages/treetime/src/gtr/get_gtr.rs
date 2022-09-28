use crate::alphabet::alphabet::AlphabetName;
use crate::gtr::gtr::GTR;
use crate::make_report;
use crate::nuc_models::jc69::{f81, hky85, jc69, k80, t92, F81Params, HKY85Params, JC69Params, K80Params, T92Params};
use clap::ArgEnum;
use eyre::{Report, WrapErr};
use smart_default::SmartDefault;
use strum_macros::Display;

#[derive(Copy, Debug, Clone, PartialEq, Eq, PartialOrd, Ord, ArgEnum, SmartDefault, Display)]
pub enum GtrModelName {
  #[default]
  JC69,
  K80,
  F81,
  HKY85,
  T92,
  Infer,
}

pub fn get_gtr(name: &GtrModelName, alphabet: &Option<AlphabetName>) -> Result<GTR, Report> {
  match name {
    GtrModelName::Infer => {
      unimplemented!("Not implemented: GTR inference is not yet implemented. Please provide the name of a concrete model with `--gtr=<model>` argument.")
    }
    GtrModelName::JC69 => jc69(JC69Params::default()),
    GtrModelName::F81 => f81(F81Params::default()),
    GtrModelName::HKY85 => hky85(HKY85Params::default()),
    GtrModelName::K80 => k80(K80Params::default()),
    GtrModelName::T92 => t92(T92Params::default()),
  }
  .wrap_err_with(|| {
    let alphabet_msg = match alphabet {
      None => "".to_owned(),
      Some(alphabet) => format!("with alphabet '{alphabet}'"),
    };
    make_report!("When creating model '{name}'{alphabet_msg}")
  })
}
