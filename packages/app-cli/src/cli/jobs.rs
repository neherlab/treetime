use clap::Args;
use serde::Serialize;
use std::thread::available_parallelism;

fn default_jobs() -> usize {
  available_parallelism().map_or(1, |n| n.get())
}

#[derive(Args, Debug, Clone)]
pub(crate) struct Jobs {
  /// Number of processing jobs. If not specified, all available CPU threads will be used.
  #[clap(
    global = true,
    display_order = 90,
    long,
    short = 'j',
    default_value_t = default_jobs(),
    hide_default_value = true
  )]
  pub jobs: usize,
}

impl Serialize for Jobs {
  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
  where
    S: serde::Serializer,
  {
    serializer.serialize_u64(self.jobs as u64)
  }
}
