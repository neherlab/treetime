use crate::representation::graph_dense::DenseGraph;
use eyre::Report;

pub fn run_optimize_dense(graph: &DenseGraph) -> Result<(), Report> {
  graph.get_edges().iter_mut().for_each(|edge| {
    let edge = edge.write_arc().payload().write_arc();
    for partition in &edge.dense_partitions {
      // TODO
    }
  });
  Ok(())
}
