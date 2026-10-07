pub mod dialect;
pub(crate) mod grammar;
pub mod model;
pub mod nexus;
pub(crate) mod nhx;
pub mod number;
pub mod read;
pub mod write;

pub use crate::dialect::{NewickAnnotations, NewickDialect, NewickStructure};
pub use crate::model::comment::{
  EdgeComment, EdgeField, LabelSide, MrBayesComment, MrBayesKind, NewickComment, NodeComment, ValueSide,
};
pub use crate::model::data::{NewickEdgeData, NewickHybrid, NewickNodeData, SupportSource};
pub use crate::model::graph::{NewickEdgeEntry, NewickGraph};
pub use crate::model::traverse::{Postorder, Preorder};
pub use crate::model::value::{NewickArray, NewickValue};
pub use crate::nexus::read::{NexusTrees, is_nexus, nexus_from_reader, nexus_from_str, nexus_trees};
pub use crate::nexus::types::{NexusCommand, NexusFile, NexusTree, NexusTreeRef, NexusWriteOptions};
pub use crate::nexus::write::{nexus_to_string, nexus_to_writer};
pub use crate::number::NumberFormat;
pub use crate::read::error::{Location, NewickError, NewickErrorKind, NewickWarning};
pub use crate::read::options::{InternalLabel, NewickReadOptions, NewickTree, ReadMode};
pub use crate::read::stream::{NewickTrees, newick_from_reader, newick_from_str, newick_trees};
pub use crate::write::conversions::{Conversion, DataKind, conversion};
pub use crate::write::newick::{newick_to_string, newick_to_writer, write_newick_trees};
pub use crate::write::options::{BranchAnnotations, NewickWriteOptions, Quoting, Spaces, SupportPlacement};

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
