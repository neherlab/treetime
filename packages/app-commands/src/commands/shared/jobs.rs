use deser::ser::Chunk;
use deser::{Atom, Error, Serialize, State};
use smart_default::SmartDefault;
use treetime_utils::cast::usize_to_u64;
use treetime_utils::init::thread_pool::available_jobs;

#[derive(Debug, Clone, SmartDefault)]
#[cfg_attr(feature = "clap", derive(clap::Args))]
pub struct JobsArgs {
  /// Number of threads of the run. Without this flag, the run uses every available CPU thread
  #[default(_code = "available_jobs()")]
  #[cfg_attr(
    feature = "clap",
    clap(
      long,
      short = 'j',
      default_value_t = available_jobs(),
      hide_default_value = true,
      help_heading = "Parallelism"
    )
  )]
  pub jobs: usize,
}

impl Serialize for JobsArgs {
  fn serialize(&self, _state: &mut State) -> Result<Chunk<'_>, Error> {
    Ok(Chunk::Atom(Atom::U64(usize_to_u64(self.jobs))))
  }
}
