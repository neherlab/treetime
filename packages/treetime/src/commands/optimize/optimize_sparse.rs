use crate::representation::graph_sparse::SparseGraph;
use eyre::Report;

pub fn run_optimize_sparse(graph: &SparseGraph) -> Result<(), Report> {
  graph.get_edges().iter_mut().for_each(|edge| {
    let mut edge = edge.write_arc().payload().write_arc();
    // for each partition, sum the parts of the likelihood messages that favor short and long branches
    for partition in &edge.sparse_partitions {
      // TODO
    }
    // balance these optimally
    edge.branch_length = Some(1.1 * edge.branch_length.unwrap_or(0.001));
  });
  Ok(())
}
