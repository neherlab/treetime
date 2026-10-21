pub mod clades;
pub mod common_ancestor;
pub mod edge;
pub mod graph;
pub mod graph_ops;
pub mod node;
pub mod node_names;
pub mod pair_by_name;
pub mod reachability;
pub mod reroot;
pub mod schedule;
pub mod schedule_scan;
pub mod topology_order;

#[cfg(test)]
mod __tests__;

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
