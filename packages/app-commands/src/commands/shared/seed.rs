use schemars::JsonSchema;
use serde::{Deserialize, Serialize};
use serde_with::skip_serializing_none;
use treetime::progress::LogSink;
use treetime::progress_info;

/// Seed of the random number generator, shared by every command that has a random step.
#[skip_serializing_none]
#[derive(Debug, Clone, Default, Serialize, Deserialize, JsonSchema)]
#[serde(default, deny_unknown_fields)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct SeedArgs {
  /// Random seed
  ///
  /// Without a seed, a run with a random step draws one and logs it, so the run can be reproduced.
  #[cfg_attr(
    feature = "clap",
    clap(long, visible_alias = "rng-seed", help_heading = "Reproducibility")
  )]
  pub seed: Option<u64>,
}

impl SeedArgs {
  pub fn resolve(&self, random_step: Option<&str>, log: &dyn LogSink) -> u64 {
    let seed = self.seed.unwrap_or_else(rand::random);
    if let Some(step) = random_step {
      progress_info!(
        log,
        "{step} is stochastic; seed {seed} (pass --seed to reproduce this run)"
      );
    }
    seed
  }
}
