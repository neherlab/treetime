use deser::{Deserialize, Serialize};
use schemars::JsonSchema;

/// Model parameter that the partitions with an inferred model share.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize, JsonSchema)]
#[deser(rename_all = "kebab-case")]
#[schemars(rename_all = "kebab-case")]
pub enum SharedParam {
  /// Exchange rates `W`, shared up to one scale per partition, so that the rate matrix of each partition has mean rate 1 under its own `pi`.
  ExchangeRates,
  /// Equilibrium frequencies `pi`.
  Frequencies,
  /// Shape alpha of the gamma rate categories.
  GammaShape,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SharedParams {
  pub exchange_rates: bool,
  pub frequencies: bool,
  pub gamma_shape: bool,
}

impl SharedParams {
  pub fn from_params(params: impl IntoIterator<Item = SharedParam>) -> Self {
    params.into_iter().fold(Self::default(), |shared, param| match param {
      SharedParam::ExchangeRates => Self {
        exchange_rates: true,
        ..shared
      },
      SharedParam::Frequencies => Self {
        frequencies: true,
        ..shared
      },
      SharedParam::GammaShape => Self {
        gamma_shape: true,
        ..shared
      },
    })
  }

  pub fn params(self) -> Vec<SharedParam> {
    [
      (SharedParam::ExchangeRates, self.exchange_rates),
      (SharedParam::Frequencies, self.frequencies),
      (SharedParam::GammaShape, self.gamma_shape),
    ]
    .into_iter()
    .filter_map(|(param, shared)| shared.then_some(param))
    .collect()
  }

  pub const fn shares_rate_matrix(self) -> bool {
    self.exchange_rates || self.frequencies
  }
}
