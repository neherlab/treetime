//! Output resolution and encoding for TreeTime application adapters.
//!
//! Owns the path resolution and format projection that turn core domain results into files: output
//! selection, path templates, compression, and the encoders for each output format. Sits above the
//! `treetime` core and below the adapters; the core does not depend on it, and adapters route their
//! output policy through it. `app-cli` is the first consumer; the deferred clients adopt it when
//! they are re-wired.

pub mod ancestral_result;
pub mod augur_node_data;
pub mod clock_result;
pub mod coalescent;
pub mod confidence;
pub mod date_comment;
pub mod mutation_comment;
pub mod optimize_result;
pub mod prune_result;
pub mod rtt;
pub mod timetree_result;

pub use date_comment::DateCommentProvider;
pub use mutation_comment::EdgeMutationCommentProvider;
pub use timetree_result::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps, TimetreeResult};

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
