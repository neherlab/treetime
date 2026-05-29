mod ancestral;
mod annotations;
mod common;
mod refine;
mod traits;

pub use crate::ancestral::{
  AugurNodeDataJsonAncestral, AugurNodeDataJsonAncestralMeta, AugurNodeDataJsonAncestralNode,
};
pub use crate::annotations::{
  AugurNodeDataJsonAnnotationEntry, AugurNodeDataJsonAnnotationSegment, AugurNodeDataJsonAnnotations,
};
pub use crate::common::{AugurNodeDataJson, AugurNodeDataJsonGeneratedBy};
pub use crate::refine::{
  AugurNodeDataJsonClock, AugurNodeDataJsonRefine, AugurNodeDataJsonRefineMeta, AugurNodeDataJsonRefineNode,
};
pub use crate::traits::{
  AugurNodeDataJsonTraitModel, AugurNodeDataJsonTraits, AugurNodeDataJsonTraitsBranches, AugurNodeDataJsonTraitsMeta,
  AugurNodeDataJsonTraitsNode,
};

#[cfg(test)]
mod __tests__;

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
