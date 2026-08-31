pub mod assign_node_names;
pub mod common_ancestor;
pub mod dependency_queue;
pub mod edge;
pub mod graph;
pub mod graph_ops;
pub mod graph_traverse;
pub mod node;
pub mod pass;
pub mod reachability;
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
