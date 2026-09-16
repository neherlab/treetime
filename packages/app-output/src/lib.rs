//! Output resolution and encoding for TreeTime application adapters.
//!
//! Owns the path resolution and format projection that turn core domain results into files: output
//! selection, path templates, compression, and the encoders for each output format. Sits above the
//! `treetime` core and below the adapters; the core does not depend on it, and adapters route their
//! output policy through it. `app-cli` is the first consumer; the deferred clients adopt it when
//! they are re-wired.

pub mod ancestral_result;
pub mod augur_node_data;
pub mod coalescent;
pub mod date_comment;
pub mod mutation_comment;
pub mod optimize_result;
pub mod prune_result;
pub mod timetree_result;

pub use date_comment::DateCommentProvider;
pub use mutation_comment::EdgeMutationCommentProvider;
pub use timetree_result::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps, TimetreeResult};
