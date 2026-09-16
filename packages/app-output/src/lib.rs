//! Output resolution and encoding for TreeTime application adapters.
//!
//! Owns the path resolution and format projection that turn core domain results into files: output
//! selection, path templates, compression, and the encoders for each output format. Sits above the
//! `treetime` core and below the adapters; the core does not depend on it, and adapters route their
//! output policy through it. `app-cli` is the first consumer; the deferred clients adopt it when
//! they are re-wired.

pub mod mutation_comment;

pub use mutation_comment::EdgeMutationCommentProvider;
