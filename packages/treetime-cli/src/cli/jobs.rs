use clap::Args;

#[derive(Args, Debug, Clone)]
pub struct Jobs {
  /// Number of processing jobs. If not specified, all available CPU threads will be used.
  #[clap(global = true, display_order = 90, long, short = 'j', default_value_t = num_cpus::get())]
  pub jobs: usize,
}
