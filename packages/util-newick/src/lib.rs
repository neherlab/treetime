pub mod annotation;
mod equality;
pub mod nexus;
mod number;
pub mod parse;
pub mod types;
mod validate;
pub mod write;

pub use crate::annotation::{write_beast_attrs, write_nhx_attrs};
pub use crate::nexus::{is_nexus, nexus_from_reader, nexus_from_string, nexus_to_string, nexus_to_writer};
pub use crate::parse::{newick_from_reader, newick_from_string};
pub use crate::types::{
  NewickEdgeData, NewickEdgeEntry, NewickGraph, NewickHybrid, NewickLabel, NewickNodeData, NewickReadOptions,
  NewickValue, NewickWriteOptions, NexusTree, NwkStyle,
};
pub use crate::write::{needs_quoting, newick_to_string, newick_to_writer, write_label};

#[cfg(test)]
mod __tests__;

#[cfg(test)]
mod tests {
  use ctor::ctor;

  #[ctor(unsafe)]
  fn init() {
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
