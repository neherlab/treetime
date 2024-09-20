use crate::representation::graph_sparse::SparseGraph;
use eyre::Report;

pub fn run_optimize_sparse(graph: &SparseGraph) -> Result<(), Report> {
  graph.get_edges().iter_mut().for_each(|edge| {
    let edge = edge.write_arc().payload().write_arc();
    for partition in &edge.sparse_partitions {
      // TODO
    }
  });
  Ok(())
}
