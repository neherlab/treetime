use deser::{Deserialize, Serialize};
use schemars::JsonSchema;
use treetime::progress::LogSink;
use treetime::progress_info;
use treetime_schema::{schema_defaults, skip_serializing_optionals};

/// Seed of the random number generator, shared by every command that has a random step.
#[derive(Debug, Clone, Default, JsonSchema, Serialize, Deserialize)]
#[deser(skip_serializing_optionals)]
#[schemars(transform = skip_serializing_optionals)]
#[schemars(default, deny_unknown_fields)]
#[schemars(transform = schema_defaults::<Self>)]
#[deser(default, deny_unknown_fields)]
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
