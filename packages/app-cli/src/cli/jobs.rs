use clap::Args;
use deser::ser::Chunk;
use deser::{Atom, Error, Serialize, State};
use treetime_utils::init::thread_pool::available_jobs;

#[derive(Args, Debug, Clone)]
pub(crate) struct Jobs {
  /// Number of processing jobs. If not specified, all available CPU threads will be used.
  #[clap(
    global = true,
    display_order = 90,
    long,
    short = 'j',
    default_value_t = available_jobs(),
    hide_default_value = true
  )]
  pub jobs: usize,
}

impl Serialize for Jobs {
  #[allow(
    clippy::as_conversions,
    reason = "count/index numeric cast is exact for the domain range"
  )]
  fn serialize(&self, _state: &mut State) -> Result<Chunk<'_>, Error> {
    Ok(Chunk::Atom(Atom::U64(self.jobs as u64)))
  }
}
