#[cfg(test)]
mod __tests__;

pub mod ancestral_result;
pub mod annotated_graph;
pub mod augur_node_data;
pub mod augur_node_data_ancestral;
pub mod augur_node_data_mugration;
pub mod augur_node_data_optimize;
pub(crate) mod auspice;
pub mod mugration_result;
pub mod mutation_filter;
pub(crate) mod nwk_comments;
pub mod optimize_result;
pub mod output_plan;
pub mod table_output;
pub mod timetree_result;
pub mod timetree_trace;
pub(crate) mod trait_profile;
pub mod tree_output;
pub(crate) mod usher_mat;

pub use timetree_result::{TimetreeEdgeOut, TimetreeNodeOut, TimetreeOutputMaps};

#[cfg(test)]
mod tests {
  use ctor::ctor;
  use treetime_utils::init::global::global_init;

  #[ctor(unsafe)]
  fn init() {
    global_init();
    rayon::ThreadPoolBuilder::new()
      .num_threads(1)
      .build_global()
      .expect("rayon global thread pool initialization failed");
  }
}
