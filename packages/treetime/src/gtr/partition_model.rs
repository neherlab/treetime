use crate::gtr::gamma::GammaRates;
use crate::gtr::gtr::GTR;
use eyre::Report;
use getset::{CopyGetters, Getters, MutGetters};
use itertools::Itertools;

#[derive(Clone, Debug, Getters, CopyGetters, MutGetters)]
pub struct PartitionModel {
  #[getset(get = "pub", get_mut = "pub(crate)")]
  gtr: GTR,
  #[getset(get = "pub")]
  gamma: GammaRates,
  #[getset(get_copy = "pub")]
  shares_gamma_shape: bool,
}

impl PartitionModel {
  pub fn new(gtr: GTR, gamma: GammaRates) -> Self {
    Self {
      gtr,
      gamma,
      shares_gamma_shape: false,
    }
  }

  pub fn without_rate_variation(gtr: GTR) -> Self {
    Self::new(gtr, GammaRates::default())
  }

  pub fn categories(&self) -> usize {
    self.gamma.categories()
  }

  #[must_use]
  pub fn with_shared_shape(self, shares_gamma_shape: bool) -> Self {
    Self {
      shares_gamma_shape,
      ..self
    }
  }

  #[must_use]
  pub fn with_gamma(self, gamma: GammaRates) -> Self {
    Self { gamma, ..self }
  }

  #[must_use]
  pub fn with_gtr(self, gtr: GTR) -> Self {
    Self { gtr, ..self }
  }

  pub fn category_gtrs(&self) -> Result<Vec<GTR>, Report> {
    self
      .gamma
      .rates()
      .iter()
      .map(|&rate| {
        let mut gtr = self.gtr.clone();
        gtr.set_mu(self.gtr.mu() * rate)?;
        Ok(gtr)
      })
      .try_collect()
  }
}

impl From<GTR> for PartitionModel {
  fn from(gtr: GTR) -> Self {
    Self::without_rate_variation(gtr)
  }
}
