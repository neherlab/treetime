pub mod nexus;
pub mod parse;
pub mod types;
pub mod write;

pub use crate::nexus::{nexus_from_reader, nexus_from_string, nexus_to_string, nexus_to_writer};
pub use crate::parse::{newick_from_reader, newick_from_string};
pub use crate::types::{
  NewickEdgeData, NewickEdgeEntry, NewickGraph, NewickHybrid, NewickNodeData, NewickValue, NewickWriteOptions,
  NexusTree, NwkStyle,
};
pub use crate::write::{
  needs_quoting, newick_to_string, newick_to_writer, write_beast_attrs, write_label, write_nhx_attrs,
};

#[cfg(test)]
mod __tests__;

#[cfg(test)]
mod tests {
  use ctor::ctor;

  #[ctor]
  fn init() {
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
