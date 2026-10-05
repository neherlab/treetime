#[cfg(test)]
mod __tests__;

pub mod annotated_graph;
pub mod augur_node_data_ancestral;
pub mod augur_node_data_refine;
pub mod augur_node_data_traits;
pub(crate) mod auspice;
pub mod mutation_filter;
pub(crate) mod nwk_comments;
pub mod output_plan;
pub mod table_output;
pub mod timetree_trace;
pub(crate) mod trait_profile;
pub mod trait_tables;
pub mod tree_output;
pub(crate) mod usher_mat;

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
