pub mod assign_node_names;
pub mod breadth_first;
pub mod common_ancestor;
pub mod edge;
pub mod find_paths;
pub mod graph;
pub mod graph_ops;
pub mod graph_traverse;
pub mod node;
pub mod reroot;
pub mod topology_order;

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
