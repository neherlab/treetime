#[cfg(test)]
mod __tests__;

pub mod ancestral_result;
pub mod ancestral_tree_output;
pub mod augur_node_data;
pub mod augur_node_data_ancestral;
pub mod augur_node_data_mugration;
pub mod augur_node_data_optimize;
pub mod clock_result;
pub mod clock_tree_output;
pub mod coalescent;
pub mod confidence;
pub mod date_comment;
pub mod discrete_trait_comment;
pub mod mugration_result;
pub mod mugration_tree_output;
pub mod mutation_comment;
pub mod optimize_result;
pub mod optimize_tree_output;
pub mod output_plan;
pub mod prune_result;
pub mod prune_tree_output;
pub mod rtt;
pub mod timetree_result;
pub mod timetree_trace;
pub mod timetree_tree_output;
pub mod tree_output;

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
