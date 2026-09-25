use eyre::Report;
use rayon::ThreadPoolBuilder;
use std::thread::available_parallelism;

pub fn available_jobs() -> usize {
  available_parallelism().map_or(1, |n| n.get())
}

pub fn init_thread_pool(jobs: usize) -> Result<(), Report> {
  let builder = ThreadPoolBuilder::new().num_threads(jobs);
  let builder = if jobs == 1 {
    builder.use_current_thread()
  } else {
    builder
  };
  builder.build_global()?;
  Ok(())
}
